function [C_opt, E_opt, Xopt, Yopt, Jc, Aopt, Bopt, Uopt, Vopt] = gridGenerator(grid_fun, Ns, Nt, Nbasis)
% gridGenerator (generalized solver only) + stochastic coordinate descent for A,B
% 1) Optimize C using Brent on objective with A=B=0.
% 2) Optimize A,B via stochastic coordinate FD updates (SGD-like).
% Objective:
%   E_total = E_NH + gammaOrth*E_orth + gammaAR*E_ar + betaAB*(||A||^2+||B||^2)
%
% Also plots optimal reference grid U,V.

%% ===== Settings =====
I = Ns; J = Nt;

% Neo-Hookean parameters
mu = 1.0;
lambda = 1.0;
epsJ = 1e-6;
kappa = 10.0;

% Spacing
hs = 1/(I-1);
ht = 1/(J-1);

% Reference denom safety in generalized solver
epsDen = 1e-14;

% Brent bracket for C
C_bracket = [0.01, 10];
C_tol = 1e-6;

% Penalty weights (tune these)
gammaOrth = 1e-2;   % orthogonality penalty strength
gammaAR   = 1e-2;   % aspect ratio penalty strength
betaAB    = 1e-4;   % regularize A,B

% SGD-like coordinate descent settings
sgd.maxIter      = 600;   % more iters; each is cheaper than full FD gradient
sgd.batchK       = 6;     % coordinates per iteration (2..12 typical)
sgd.fdStep       = 1e-4;  % finite difference step
sgd.lr0          = 0.08;  % learning rate start (try 0.03..0.2)
sgd.decay        = 1e-3;  % lr = lr0/(1+decay*it)
sgd.displayEvery = 25;

%% ===== Boundary setup =====
[Xc, Yc] = grid_fun(I-1, J-1);

BND.X = nan(J,I); BND.Y = nan(J,I);
BND.X(1,:) = Xc(1,:);    BND.Y(1,:) = Yc(1,:);
BND.X(J,:) = Xc(end,:);  BND.Y(J,:) = Yc(end,:);
BND.X(:,1) = Xc(:,1);    BND.Y(:,1) = Yc(:,1);
BND.X(:,I) = Xc(:,end);  BND.Y(:,I) = Yc(:,end);

% corners
BND.X(1,1)=Xc(1,1);         BND.Y(1,1)=Yc(1,1);
BND.X(1,I)=Xc(1,end);       BND.Y(1,I)=Yc(1,end);
BND.X(J,1)=Xc(end,1);       BND.Y(J,1)=Yc(end,1);
BND.X(J,I)=Xc(end,end);     BND.Y(J,I)=Yc(end,end);

%% ===== Bubble basis (m+n ordering): sin(pi*m*s)sin(pi*n*t) =====
[sgrid, tgrid] = meshgrid(linspace(0,1,I), linspace(0,1,J));
modes = makeModes_mplusn(Nbasis);
basisFns = cell(Nbasis,1);
for k = 1:Nbasis
    m = modes(k,1); n = modes(k,2);
    basisFns{k} = sin(pi*m*sgrid).*sin(pi*n*tgrid);
end

%% ===== Stage 1: Optimize C using generalized solver with A=B=0 =====
A0 = zeros(Nbasis,1);
B0 = zeros(Nbasis,1);

fC = @(C) objectiveGeneralized(C, A0, B0);
[C_opt, ~] = brentMin(fC, C_bracket(1), C_bracket(2), C_tol);
fprintf('C_opt = %.8f (optimized using generalized solver, A=B=0)\n', C_opt);

%% ===== Stage 2: Optimize A,B using stochastic coordinate FD updates =====
if Nbasis == 0
    Aopt = zeros(0,1); Bopt = zeros(0,1);
else
    z0 = zeros(2*Nbasis,1); % z=[A;B]
    fAB = @(z) objectiveGeneralized(C_opt, z(1:Nbasis), z(Nbasis+1:end));
    [zopt, E_opt] = sgdCoordFD(fAB, z0, sgd);
    Aopt = zopt(1:Nbasis);
    Bopt = zopt(Nbasis+1:end);
end

% Final grid + final energy (computed at final A,B)
[Xopt, Yopt] = solveGeneralizedLengthGrid(C_opt, BND, Aopt, Bopt, basisFns, epsDen);
E_opt = objectiveGeneralized(C_opt, Aopt, Bopt);

fprintf('Final E_opt (NH + penalties) = %.6e\n', E_opt);

%% ===== Jacobian at cell centers =====
[xs, ys, xt, yt] = metricsCentered(Xopt, Yopt, hs, ht);
Jc = xs.*yt - xt.*ys;
fprintf('min(Jc) = %.6e\n', min(Jc(:)));

%% ===== Build optimal reference grid =====
[Uopt, Vopt] = buildReferenceGrid(C_opt, Aopt, Bopt, basisFns, I, J);

%% ===== Plot final physical grid =====
figure; hold on; axis equal; box on;
for j=1:J, plot(Xopt(j,:), Yopt(j,:), '-k'); end
for i=1:I, plot(Xopt(:,i), Yopt(:,i), '-k'); end
title(sprintf('Physical grid (C*=%.4f, Nbasis=%d)', C_opt, Nbasis));
xlabel('X'); ylabel('Y');

%% ===== Plot Jacobian heatmap =====
figure; imagesc(Jc); axis image; colorbar;
title('Jacobian det J at cell centers'); xlabel('i'); ylabel('j');

%% ===== Plot optimal reference grid (u,v) =====
figure; hold on; axis equal; box on;
for j=1:J, plot(Uopt(j,:), Vopt(j,:), '-b'); end
for i=1:I, plot(Uopt(:,i), Vopt(:,i), '-b'); end
title(sprintf('Reference grid (u,v), C*=%.4f, Nbasis=%d', C_opt, Nbasis));
xlabel('u'); ylabel('v');

%% ===== Overlay reference vs physical =====
figure; hold on; axis equal; box on;
for j=1:J, plot(Uopt(j,:), Vopt(j,:), 'b-'); end
for i=1:I, plot(Uopt(:,i), Vopt(:,i), 'b-'); end
for j=1:J, plot(Xopt(j,:), Yopt(j,:), 'k-'); end
for i=1:I, plot(Xopt(:,i), Yopt(:,i), 'k-'); end
legend('ref rows','ref cols','phys rows','phys cols');
title('Reference (blue) vs Physical (black)');
xlabel('x/u'); ylabel('y/v');

%% ===== nested objective =====
    function E = objectiveGeneralized(C, A, B)
        [X, Y] = solveGeneralizedLengthGrid(C, BND, A, B, basisFns, epsDen);

        % Base NH energy
        E_NH = neoHookeanEnergy(X, Y, hs, ht, mu, lambda, epsJ, kappa);

        % Penalties
        E_orth = orthogonalityPenalty(X, Y, hs, ht);
        E_ar   = aspectRatioPenalty(X, Y, hs, ht);

        % Total (plus mild AB regularization)
        E = E_NH + gammaOrth*E_orth + gammaAR*E_ar + betaAB*(sum(A.^2)+sum(B.^2));
    end
end


%% =======================================================================
%% Build reference grid u,v on parameter grid
function [U, V] = buildReferenceGrid(C, A, B, basisFns, I, J)
    [sgrid, tgrid] = meshgrid(linspace(0,1,I), linspace(0,1,J));
    U = sgrid;
    V = C*tgrid;
    for k = 1:numel(A)
        if isempty(basisFns{k}), continue; end
        U = U + A(k)*basisFns{k};
        V = V + B(k)*basisFns{k};
    end
end


%% =======================================================================
%% Generalized-length solver (denominator squared), SAME operator for X and Y
function [X, Y] = solveGeneralizedLengthGrid(C, BND, A, B, basisFns, epsDen)
    [J,I] = size(BND.X);
    nI = I-2; nJ = J-2;
    N  = nI*nJ;

    [sgrid, tgrid] = meshgrid(linspace(0,1,I), linspace(0,1,J));
    u = sgrid;
    v = C*tgrid;

    for k=1:numel(A)
        if isempty(basisFns{k}), continue; end
        u = u + A(k)*basisFns{k};
        v = v + B(k)*basisFns{k};
    end

    du_h = u(:,2:end) - u(:,1:end-1);    % J x (I-1)
    dv_h = v(:,2:end) - v(:,1:end-1);
    du_v = u(2:end,:) - u(1:end-1,:);    % (J-1) x I
    dv_v = v(2:end,:) - v(1:end-1,:);

    % Smooth floor: den = den + epsDen (avoids kinks from max())
    den_h = (du_h.^2 + dv_h.^2) + epsDen;
    den_v = (du_v.^2 + dv_v.^2) + epsDen;

    Wh = 1 ./ den_h;
    Wv = 1 ./ den_v;

    idx = @(ii,jj) (jj-1)*nI + ii;

    ii = zeros(5*N,1);
    jj = zeros(5*N,1);
    ss = zeros(5*N,1);
    bX = zeros(N,1);
    bY = zeros(N,1);
    ptr = 0;

    for jj_int = 1:nJ
        jg = jj_int + 1;
        for ii_int = 1:nI
            ig  = ii_int + 1;
            row = idx(ii_int, jj_int);

            Ww = Wh(jg, ig-1);
            We = Wh(jg, ig);
            Ws = Wv(jg-1, ig);
            Wn = Wv(jg, ig);

            diagv = Ww + We + Ws + Wn;

            ptr=ptr+1; ii(ptr)=row; jj(ptr)=row; ss(ptr)=diagv;

            % west
            if ii_int > 1
                ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int-1,jj_int); ss(ptr)=-Ww;
            else
                bX(row) = bX(row) + Ww * BND.X(jg,1);
                bY(row) = bY(row) + Ww * BND.Y(jg,1);
            end

            % east
            if ii_int < nI
                ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int+1,jj_int); ss(ptr)=-We;
            else
                bX(row) = bX(row) + We * BND.X(jg,I);
                bY(row) = bY(row) + We * BND.Y(jg,I);
            end

            % south
            if jj_int > 1
                ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int,jj_int-1); ss(ptr)=-Ws;
            else
                bX(row) = bX(row) + Ws * BND.X(1,ig);
                bY(row) = bY(row) + Ws * BND.Y(1,ig);
            end

            % north
            if jj_int < nJ
                ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int,jj_int+1); ss(ptr)=-Wn;
            else
                bX(row) = bX(row) + Wn * BND.X(J,ig);
                bY(row) = bY(row) + Wn * BND.Y(J,ig);
            end
        end
    end

    Aop = sparse(ii(1:ptr), jj(1:ptr), ss(1:ptr), N, N);

    Xvec = Aop \ bX;
    Yvec = Aop \ bY;

    X = BND.X; Y = BND.Y;
    X(2:end-1,2:end-1) = reshape(Xvec, [nI,nJ])';
    Y(2:end-1,2:end-1) = reshape(Yvec, [nI,nJ])';
end


%% =======================================================================
%% Neo-Hookean energy (2D) with barrier
function E = neoHookeanEnergy(X, Y, hs, ht, mu, lambda, epsJ, kappa)
    cellArea = hs*ht;

    [xs, ys, xt, yt] = metricsCentered(X, Y, hs, ht);
    Jc = xs.*yt - xt.*ys;

    Jsafe = max(Jc, epsJ);
    logJ  = log(Jsafe);

    d = 2;
    trC = xs.^2 + ys.^2 + xt.^2 + yt.^2;
    W = 0.5*mu*(trC - d) - mu*logJ + 0.5*lambda*(logJ.^2);

    r = max(0, epsJ - Jc);
    B = kappa*(r.^2);

    E = cellArea * sum(W(:) + B(:));
end


%% =======================================================================
%% Orthogonality penalty: ∫ (Xs·Xt)^2
function Eorth = orthogonalityPenalty(X, Y, hs, ht)
    [xs, ys, xt, yt] = metricsCentered(X, Y, hs, ht);
    orth = xs.*xt + ys.*yt;
    Eorth = (hs*ht) * sum(orth(:).^2);
end


%% =======================================================================
%% Aspect ratio penalty: ∫ log(||Xs||^2 / ||Xt||^2)^2
function Ear = aspectRatioPenalty(X, Y, hs, ht)
    [xs, ys, xt, yt] = metricsCentered(X, Y, hs, ht);
    a = xs.^2 + ys.^2;
    b = xt.^2 + yt.^2;
    epsr = 1e-12;
    r = log((a+epsr)./(b+epsr));
    Ear = (hs*ht) * sum(r(:).^2);
end


%% =======================================================================
%% Centered metrics at cell centers (J-2 x I-2)
function [xs, ys, xt, yt] = metricsCentered(X, Y, hs, ht)
    x_s = (X(:,3:end)-X(:,1:end-2))/(2*hs);
    y_s = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
    x_t = (X(3:end,:)-X(1:end-2,:))/(2*ht);
    y_t = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);

    xs = x_s(2:end-1,:);
    ys = y_s(2:end-1,:);
    xt = x_t(:,2:end-1);
    yt = y_t(:,2:end-1);
end


%% =======================================================================
%% Modes ordered by m+n (low frequency first)
function modes = makeModes_mplusn(N)
    if N <= 0
        modes = zeros(0,2);
        return;
    end
    modes = zeros(N,2);
    k = 0;
    nsum = 2;
    while k < N
        for m = 1:(nsum-1)
            n = nsum - m;
            if n >= 1
                k = k + 1;
                modes(k,:) = [m n];
                if k == N, break; end
            end
        end
        nsum = nsum + 1;
    end
end


%% =======================================================================
%% Stochastic coordinate descent with FD gradients (SGD-like)
function [z, fz] = sgdCoordFD(f, z0, opts)
    if ~isfield(opts,'maxIter'), opts.maxIter = 500; end
    if ~isfield(opts,'batchK'),  opts.batchK  = 6; end
    if ~isfield(opts,'fdStep'),  opts.fdStep  = 1e-4; end
    if ~isfield(opts,'lr0'),     opts.lr0     = 0.1; end
    if ~isfield(opts,'decay'),   opts.decay   = 1e-3; end
    if ~isfield(opts,'displayEvery'), opts.displayEvery = 25; end

    z = z0(:);
    fz = f(z);

    n = numel(z);
    K = min(opts.batchK, n);

    for it = 1:opts.maxIter
        lr = opts.lr0 / (1 + opts.decay*it);

        idxs = randperm(n, K);
        ghat = zeros(n,1);

        for kk = 1:K
            i = idxs(kk);
            hi = opts.fdStep * max(1, abs(z(i)));
            zp = z; zm = z;
            zp(i) = zp(i) + hi;
            zm(i) = zm(i) - hi;
            ghat(i) = (f(zp) - f(zm)) / (2*hi);
        end

        z_new = z - lr*ghat;
        f_new = f(z_new);

        if f_new <= fz
            z = z_new; fz = f_new;
        else
            % one-step shrink safeguard
            lr2 = 0.5*lr;
            z_new = z - lr2*ghat;
            f_new = f(z_new);
            if f_new <= fz
                z = z_new; fz = f_new;
            end
        end

        if mod(it, opts.displayEvery) == 0
            fprintf('SCD iter %4d: f=%.6e  lr=%.3e  K=%d\n', it, fz, lr, K);
        end
    end
end


%% =======================================================================
%% Brent minimizer for scalar f on [a,b]
function [xopt, fopt] = brentMin(f, a, b, tol)
    if nargin < 4, tol = 1e-6; end
    phi  = (3 - sqrt(5))/2;
    eps0 = 1e-12;

    x = a + phi*(b-a);
    w = x; v = x;
    fx = f(x); fw = fx; fv = fx;
    d = 0; e = 0;

    while (b-a) > tol
        m = 0.5*(a+b);
        tol1 = tol*abs(x) + eps0;
        tol2 = 2*tol1;

        if abs(x-m) <= (tol2 - 0.5*(b-a))
            break;
        end

        p = 0; q = 0; r = 0;

        if abs(x-w) > eps0 && abs(x-v) > eps0 && abs(w-v) > eps0
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q - (x-w)*r;
            q = 2*(q-r);
            if q > 0, p = -p; end
            q = abs(q);

            if abs(p) < abs(0.5*q*e) && p > q*(a-x) && p < q*(b-x)
                d = p/q;
                u = x + d;
                if (u-a) < tol2 || (b-u) < tol2
                    d = sign(m-x)*tol1;
                end
            else
                if x < m, e = b-x; else, e = a-x; end
                d = phi*e;
            end
        else
            if x < m, e = b-x; else, e = a-x; end
            d = phi*e;
        end

        if abs(d) >= tol1
            u = x + d;
        else
            u = x + sign(d)*tol1;
        end

        fu = f(u);

        if fu <= fx
            if u < x, b = x; else, a = x; end
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        else
            if u < x, a = u; else, b = u; end
            if fu <= fw || w == x
                v = w; fv = fw;
                w = u; fw = fu;
            elseif fu <= fv || v == x || v == w
                v = u; fv = fu;
            end
        end
    end

    xopt = x;
    fopt = fx;
end
