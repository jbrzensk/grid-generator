function [C_opt, E_opt, Xopt, Yopt, Jc, pUopt, pVopt, Uopt, Vopt, info] = gridGenerator(grid_fun, Ns, Nt, Nsteps)
% gridGenerator (RECLUSTER LOOP VERSION - DUAL-CELL CONSISTENT)
% Cheap robust reclustering on the reference rectangle using a Jacobian-based monitor.
% Nsteps = number of reclustering iterations.
%
% This version avoids padarray and uses widths consistent with the monitor layout:
%   - if w is (J-2)x(I-2), use dual-cell widths
%   - if w is (J-1)x(I-1), use full strip widths
%   - if w is JxI, first average to cell values

%% ---- initialize outputs ----
C_opt = NaN; E_opt = NaN;
Xopt = []; Yopt = []; Jc = [];
pUopt = []; pVopt = [];
Uopt = []; Vopt = [];
info = struct();

%% ---------------- settings ----------------
I = Ns; J = Nt;
if nargin < 4 || isempty(Nsteps)
    Nsteps = 25;
end
Nsteps = floor(Nsteps);

derivOrder = 4;

thetaW = 0.2;
mu = 1.0; lambda = 1.0;
epsJ = 1e-6; kappa = 10.0;
epsJw = 1e-10;

gammaOrth = 1e-2;
gammaAR   = 1e-2;
betaAB    = 1e-6;

filterMode = 'minq';
qminSJ   = 0.25;

epsDen = 1e-14;

hs = 1/(I-1);
ht = 1/(J-1);

useBrent = true;
C_fixed = 0.4269;
C_bracket = [0.01, 10];
C_tol = 1e-6;

alphaMon = 0.5;
epsMon   = 1e-12;
doU = true;
doV = true;
nSmooth = 1;

%% ---------------- boundary setup ----------------
[Xc, Yc] = grid_fun(I-1, J-1);

if ~isequal(size(Xc), [J I]) && isequal(size(Xc), [I J])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.'; Yc = Yc.';
end
if ~isequal(size(Xc), [J I])
    error('grid_fun must return size (Nt x Ns) = (%d x %d). Got %s.', J, I, mat2str(size(Xc)));
end

BND.X = nan(J,I); BND.Y = nan(J,I);
BND.X(1,:) = Xc(1,:);    BND.Y(1,:) = Yc(1,:);
BND.X(J,:) = Xc(end,:);  BND.Y(J,:) = Yc(end,:);
BND.X(:,1) = Xc(:,1);    BND.Y(:,1) = Yc(:,1);
BND.X(:,I) = Xc(:,end);  BND.Y(:,I) = Yc(:,end);

BND.X(1,1)=Xc(1,1);         BND.Y(1,1)=Yc(1,1);
BND.X(1,I)=Xc(1,end);       BND.Y(1,I)=Yc(1,end);
BND.X(J,1)=Xc(end,1);       BND.Y(J,1)=Yc(end,1);
BND.X(J,I)=Xc(end,end);     BND.Y(J,I)=Yc(end,end);

%% ---------------- Stage 1: choose C ----------------
pU0 = zeros(I-1,1);
pV0 = zeros(J-1,1);

if useBrent
    fC = @(C) objectiveScalar(C, pU0, pV0);
    [C_opt, ~] = brentMin(fC, C_bracket(1), C_bracket(2), C_tol);
else
    C_opt = C_fixed;
end
fprintf('C_opt = %.8f\n', C_opt);

%% ---------------- Stage 2: reclustering loop ----------------
pU = pU0;
pV = pV0;

Ehist = nan(max(Nsteps,1),1);
CVhist = nan(max(Nsteps,1),1);
Vloghist = nan(max(Nsteps,1),1);

if Nsteps <= 0
    [Uopt, Vopt] = buildReferenceGridSoftmax(C_opt, pU, pV, I, J);
    [Xopt, Yopt] = solveAnisoLaplaceFromUV(Uopt, Vopt, BND, epsDen, hs, ht, derivOrder);
    [xs, ys, xt, yt] = metricsHighOrder(Xopt, Yopt, hs, ht, derivOrder);
    Jc = xs.*yt - xt.*ys;

    pUopt = pU; pVopt = pV;
    E_opt = NaN;
    info.Ehist = [];
    info.CVhist = [];
    info.Vloghist = [];
else
    for k = 1:Nsteps
        [U, V] = buildReferenceGridSoftmax(C_opt, pU, pV, I, J);
        [X, Y] = solveAnisoLaplaceFromUV(U, V, BND, epsDen, hs, ht, derivOrder);

        ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder);
        EW  = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw);
        Ebase = (1-thetaW)*ENH + thetaW*EW;

        Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder);
        Ear   = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder);
        Ereg  = betaAB*(sum(pU(:).^2)+sum(pV(:).^2));
        Ehist(k) = Ebase + gammaOrth*Eorth + gammaAR*Ear + Ereg;

        [xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
        Jc_now = xs.*yt - xt.*ys;

        muJ = mean(Jc_now(:));
        sigJ = std(Jc_now(:));
        CVhist(k) = sigJ / max(1e-30, abs(muJ));
        Vloghist(k) = var(log(abs(Jc_now(:)) + 1e-30));

        w = (abs(Jc_now) + epsMon).^alphaMon;

        [pU, pV] = reclusterSoftmaxFromMonitor(U, V, w, C_opt, doU, doV, nSmooth);

        if mod(k, max(1, round(Nsteps/10))) == 0 || k == 1
            fprintf('recluster %4d/%d: E=%.6e  CV(J)=%.3e  Var(log|J|)=%.3e\n', ...
                k, Nsteps, Ehist(k), CVhist(k), Vloghist(k));
        end
    end

    info.Ehist = Ehist(1:Nsteps);
    info.CVhist = CVhist(1:Nsteps);
    info.Vloghist = Vloghist(1:Nsteps);

    E_opt = Ehist(Nsteps);
    pUopt = pU;
    pVopt = pV;

    [Uopt, Vopt] = buildReferenceGridSoftmax(C_opt, pUopt, pVopt, I, J);
    [Xopt, Yopt] = solveAnisoLaplaceFromUV(Uopt, Vopt, BND, epsDen, hs, ht, derivOrder);

    [xs, ys, xt, yt] = metricsHighOrder(Xopt, Yopt, hs, ht, derivOrder);
    Jc = xs.*yt - xt.*ys;
end

[minq, ESJ] = scaledJacobianStatsHO(Xopt, Yopt, hs, ht, derivOrder, qminSJ);
if strcmpi(filterMode,'minq')
    g_opt = minq;
else
    g_opt = ESJ;
end
info.minq = minq;
info.ESJ  = ESJ;
info.g_opt = g_opt;

fprintf('Final: E=%.6e, minq=%.6e, ESJ=%.6e\n', E_opt, minq, ESJ);

%% ---------------- plots ----------------
figure; hold on; axis equal; box on;
for jj=1:J, plot(Xopt(jj,:), Yopt(jj,:), '-k'); end
for ii=1:I, plot(Xopt(:,ii), Yopt(:,ii), '-k'); end
title(sprintf('Physical grid (recluster loop), C=%.4f', C_opt));
xlabel('X'); ylabel('Y');

figure; hold on; axis equal; box on;
for jj=1:J, plot(Uopt(jj,:), Vopt(jj,:), '-b'); end
for ii=1:I, plot(Uopt(:,ii), Vopt(:,ii), '-b'); end
title('Reference grid (reclustered U,V)');
xlabel('u'); ylabel('v');

figure; imagesc(Jc); axis image; colorbar;
title('det(dX/dsdt) at cell centers'); xlabel('i'); ylabel('j');

if Nsteps > 0
    figure; plot(1:Nsteps, info.CVhist, '-o'); grid on;
    xlabel('recluster iter'); ylabel('CV(detJ)');

    figure; plot(1:Nsteps, info.Vloghist, '-o'); grid on;
    xlabel('recluster iter'); ylabel('Var(log|detJ|)');
end

%% ---------------- nested objective helpers ----------------
    function E = objectiveScalar(C, pUloc, pVloc)
        [E, ~] = objectiveWithSliver(C, pUloc, pVloc); %#ok<ASGLU>
    end

    function [E, g] = objectiveWithSliver(C, pUloc, pVloc)
        [U, V] = buildReferenceGridSoftmax(C, pUloc, pVloc, I, J);
        [X, Y] = solveAnisoLaplaceFromUV(U, V, BND, epsDen, hs, ht, derivOrder);

        ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder);
        EW  = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw);
        Ebase = (1-thetaW)*ENH + thetaW*EW;

        Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder);
        Ear   = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder);
        Ereg  = betaAB*(sum(pUloc(:).^2)+sum(pVloc(:).^2));

        E = Ebase + gammaOrth*Eorth + gammaAR*Ear + Ereg;

        [minq2, ESJ2] = scaledJacobianStatsHO(X, Y, hs, ht, derivOrder, qminSJ);
        if strcmpi(filterMode,'minq')
            g = minq2;
        else
            g = ESJ2;
        end
    end
end

%% =======================================================================
function [pUnew, pVnew] = reclusterSoftmaxFromMonitor(Ugrid, Vgrid, w, C, doU, doV, nSmooth)
% Reclustering using the monitor on the SAME layout it is computed on.
% Supports:
%   w size = (J-2)x(I-2)  [preferred: from metricsHighOrder Jacobian]
%   w size = (J-1)x(I-1)
%   w size = JxI          [node-based; converted to cells]
%
% No padarray, no toolbox dependencies.

[J,I] = size(Ugrid);

U = Ugrid(1,:).';
V = Vgrid(:,1);

du = diff(U);       % (I-1)x1
dv = diff(V);       % (J-1)x1

sw = size(w);

% normalize orientation first
if isequal(sw,[I-2,J-2]) || isequal(sw,[I-1,J-1]) || isequal(sw,[I,J])
    w = w.';
    sw = size(w);
end

% convert node-based monitor to cell-based if needed
if isequal(sw,[J,I])
    w = 0.25*(w(1:end-1,1:end-1) + w(2:end,1:end-1) + ...
              w(1:end-1,2:end)   + w(2:end,2:end));
    sw = size(w); % now (J-1)x(I-1)
end

% choose widths consistent with w layout
if isequal(sw,[J-2,I-2])
    du_eff = 0.5*(du(1:end-1) + du(2:end));   % (I-2)x1
    dv_eff = 0.5*(dv(1:end-1) + dv(2:end));   % (J-2)x1
elseif isequal(sw,[J-1,I-1])
    du_eff = du;                               % (I-1)x1
    dv_eff = dv;                               % (J-1)x1
else
    error('Monitor w has unsupported size %s. Expected [%d %d], [%d %d], or [%d %d].', ...
        mat2str(sw), J-2, I-2, J-1, I-1, J, I);
end

if nSmooth > 0
    for s = 1:nSmooth
        w = smooth2D(w);
    end
end

% ---------------- V recluster ----------------
if doV
    A = w * du_eff;          % row masses
    Mv = A .* dv_eff;        % weighted strip mass

    Stot = sum(Mv);
    if ~isfinite(Stot) || Stot <= 0, Stot = 1; end

    if isequal(sw,[J-2,I-2])
        Vnodes = [V(1); V(2:end-1); V(end)];
        Snodes = [0; cumsum(Mv); Stot];
    else
        Vnodes = V;
        Snodes = [0; cumsum(Mv)];
    end

    Vnew = invertCumulative(Vnodes, Snodes, linspace(0, Stot, J).');
    Vnew(1) = 0; Vnew(end) = C;
    Vnew = enforceMonotone(Vnew);

    dvNew = diff(Vnew);
    dvNew = max(dvNew, 1e-15);
    pVnew = log(dvNew);
    pVnew = pVnew - mean(pVnew);
else
    pVnew = log(diff(V));
    pVnew = pVnew - mean(pVnew);
end

% ---------------- U recluster ----------------
if doU
    B = w.' * dv_eff;        % column masses
    Mu = B .* du_eff;

    Ttot = sum(Mu);
    if ~isfinite(Ttot) || Ttot <= 0, Ttot = 1; end

    if isequal(sw,[J-2,I-2])
        Unodes = [U(1); U(2:end-1); U(end)];
        Tnodes = [0; cumsum(Mu); Ttot];
    else
        Unodes = U;
        Tnodes = [0; cumsum(Mu)];
    end

    Unew = invertCumulative(Unodes, Tnodes, linspace(0, Ttot, I).');
    Unew(1) = 0; Unew(end) = 1;
    Unew = enforceMonotone(Unew);

    duNew = diff(Unew);
    duNew = max(duNew, 1e-15);
    pUnew = log(duNew);
    pUnew = pUnew - mean(pUnew);
else
    pUnew = log(diff(U));
    pUnew = pUnew - mean(pUnew);
end
end

function W = smooth2D(W)
W = W(:,:);
if min(size(W)) < 3, return; end
W2 = W;
W2(2:end-1,2:end-1) = 0.5*W(2:end-1,2:end-1) + ...
    0.125*(W(1:end-2,2:end-1) + W(3:end,2:end-1) + ...
           W(2:end-1,1:end-2) + W(2:end-1,3:end));
W = W2;
end

function x = enforceMonotone(x)
x = x(:);
for k = 2:numel(x)
    if x(k) < x(k-1)
        x(k) = x(k-1);
    end
end
end

function xnew = invertCumulative(xnodes, S, Squery)
xnodes = xnodes(:);
S = S(:);
Squery = Squery(:);

for k = 2:numel(S)
    if S(k) < S(k-1)
        S(k) = S(k-1);
    end
end

keep = [true; diff(S) > 0];
if sum(keep) < 2
    xnew = linspace(xnodes(1), xnodes(end), numel(Squery)).';
    return;
end

xnew = interp1(S(keep), xnodes(keep), Squery, 'linear', 'extrap');
end

%% =======================================================================
function [Ugrid, Vgrid] = buildReferenceGridSoftmax(C, pU, pV, I, J)
pU = pU(:); pV = pV(:);
dU = softmax_safe(pU);
dV = softmax_safe(pV);
dU = dU / sum(dU);
dV = dV / sum(dV) * C;
U = [0; cumsum(dU)];
V = [0; cumsum(dV)];
Ugrid = repmat(U(:).', J, 1);
Vgrid = repmat(V(:),    1, I);
end

function y = softmax_safe(x)
x = x(:);
if exist('softmax','file') == 2
    y = softmax(x);
else
    x = x - max(x);
    ex = exp(x);
    y = ex / sum(ex);
end
end

%% =======================================================================
function [minq, ESJ] = scaledJacobianStatsHO(X, Y, hs, ht, derivOrder, qmin)
if nargin<6, qmin=0.25; end
epsn = 1e-14;
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
Jc = xs.*yt - xt.*ys;
na = sqrt(xs.^2 + ys.^2 + epsn);
nb = sqrt(xt.^2 + yt.^2 + epsn);
q = Jc ./ (na.*nb + epsn);
minq = min(q(:));
r = max(0, qmin - q);
ESJ = (hs*ht) * sum(r(:).^2);
end

%% =======================================================================
function [X, Y] = solveAnisoLaplaceFromUV(Ugrid, Vgrid, BND, epsDen, hs, ht, derivOrder)
[J,I] = size(BND.X);
nI = I-2; nJ = J-2;
N  = nI*nJ;

u = Ugrid; v = Vgrid;

[usN, vsN, utN, vtN] = metricsHighOrder_nodes(u, v, hs, ht, derivOrder);

alphaN = 1 ./ (usN.^2 + vsN.^2 + epsDen);
betaN  = 1 ./ (utN.^2 + vtN.^2 + epsDen);

alpha_e = 0.5*(alphaN(:,1:I-1) + alphaN(:,2:I));
beta_n  = 0.5*(betaN(1:J-1,:) + betaN(2:J,:));

As = alpha_e / (hs^2);
At = beta_n  / (ht^2);

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

        Ww = As(jg, ig-1);
        We = As(jg, ig);
        Ws = At(jg-1, ig);
        Wn = At(jg, ig);

        diagv = Ww + We + Ws + Wn;

        ptr=ptr+1; ii(ptr)=row; jj(ptr)=row; ss(ptr)=diagv;

        if ii_int > 1
            ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int-1,jj_int); ss(ptr)=-Ww;
        else
            bX(row) = bX(row) + Ww * BND.X(jg,1);
            bY(row) = bY(row) + Ww * BND.Y(jg,1);
        end

        if ii_int < nI
            ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int+1,jj_int); ss(ptr)=-We;
        else
            bX(row) = bX(row) + We * BND.X(jg,I);
            bY(row) = bY(row) + We * BND.Y(jg,I);
        end

        if jj_int > 1
            ptr=ptr+1; ii(ptr)=row; jj(ptr)=idx(ii_int,jj_int-1); ss(ptr)=-Ws;
        else
            bX(row) = bX(row) + Ws * BND.X(1,ig);
            bY(row) = bY(row) + Ws * BND.Y(1,ig);
        end

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
function ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder)
cellArea = hs*ht;
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
Jc = xs.*yt - xt.*ys;
Jsafe = max(Jc, epsJ);
logJ  = log(Jsafe);
d = 2;
trC = xs.^2 + ys.^2 + xt.^2 + yt.^2;
W = 0.5*mu*(trC - d) - mu*logJ + 0.5*lambda*(logJ.^2);
r = max(0, epsJ - Jc);
B = kappa*(r.^2);
ENH = cellArea * sum(W(:) + B(:));
end

function EW = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw)
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
a2 = xs.^2 + ys.^2;
b2 = xt.^2 + yt.^2;
Jc = xs.*yt - ys.*xt;
Jsafe = max(Jc, epsJw);
tmp = (a2 + b2) ./ Jsafe;
EW = (hs*ht) * sum(tmp(:));
end

function Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder)
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
orth = xs.*xt + ys.*yt;
Eorth = (hs*ht) * sum(orth(:).^2);
end

function Ear = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder)
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
a = xs.^2 + ys.^2;
b = xt.^2 + yt.^2;
epsr = 1e-12;
r = log((a+epsr)./(b+epsr));
Ear = (hs*ht) * sum(r(:).^2);
end

%% =======================================================================
function [xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder)
Xs = diff1_highorder_dim2(X, hs, derivOrder);
Ys = diff1_highorder_dim2(Y, hs, derivOrder);
Xt = diff1_highorder_dim1(X, ht, derivOrder);
Yt = diff1_highorder_dim1(Y, ht, derivOrder);
xs = Xs(2:end-1,2:end-1);
ys = Ys(2:end-1,2:end-1);
xt = Xt(2:end-1,2:end-1);
yt = Yt(2:end-1,2:end-1);
end

function [xs, ys, xt, yt] = metricsHighOrder_nodes(X, Y, hs, ht, derivOrder)
xs = diff1_highorder_dim2(X, hs, derivOrder);
ys = diff1_highorder_dim2(Y, hs, derivOrder);
xt = diff1_highorder_dim1(X, ht, derivOrder);
yt = diff1_highorder_dim1(Y, ht, derivOrder);
end

function Df = diff1_highorder_dim2(F, h, p)
if mod(p,2)~=0 || p<2, error('derivOrder must be even >=2'); end
[J,I] = size(F);
r = p/2; m = 2*r + 1;
Df = zeros(J,I);
x = (-r:r);
wc = fdweights(0, x, 1);

for i = (r+1):(I-r)
    Df(:,i) = (F(:, i-r:i+r) * wc(:)) / h;
end
for i = 1:r
    xL = (0:m-1) - (i-1);
    wL = fdweights(0, xL, 1);
    Df(:,i) = (F(:, 1:m) * wL(:)) / h;
end
for i = (I-r+1):I
    xR = (-(m-1):0) - (i-I);
    wR = fdweights(0, xR, 1);
    Df(:,i) = (F(:, I-m+1:I) * wR(:)) / h;
end
end

function Df = diff1_highorder_dim1(F, h, p)
if mod(p,2)~=0 || p<2, error('derivOrder must be even >=2'); end
[J,I] = size(F);
r = p/2; m = 2*r + 1;
Df = zeros(J,I);
x = (-r:r);
wc = fdweights(0, x, 1);

for j = (r+1):(J-r)
    Df(j,:) = (wc(:).' * F(j-r:j+r, :)) / h;
end
for j = 1:r
    xB = (0:m-1) - (j-1);
    wB = fdweights(0, xB, 1);
    Df(j,:) = (wB(:).' * F(1:m, :)) / h;
end
for j = (J-r+1):J
    xT = (-(m-1):0) - (j-J);
    wT = fdweights(0, xT, 1);
    Df(j,:) = (wT(:).' * F(J-m+1:J, :)) / h;
end
end

function w = fdweights(x0, x, m)
n = numel(x);
c = zeros(n, m+1);
c1 = 1;
c4 = x(1) - x0;
c(1,1) = 1;
for i = 2:n
    mn = min(i, m+1);
    c2 = 1;
    c5 = c4;
    c4 = x(i) - x0;
    for j = 1:i-1
        c3 = x(i) - x(j);
        c2 = c2 * c3;
        if j == i-1
            for k = mn:-1:2
                c(i,k) = (c1*((k-1)*c(i-1,k-1) - c5*c(i-1,k))) / c2;
            end
            c(i,1) = (-c1*c5*c(i-1,1)) / c2;
        end
        for k = mn:-1:2
            c(j,k) = ((c4*c(j,k)) - (k-1)*c(j,k-1)) / c3;
        end
        c(j,1) = (c4*c(j,1)) / c3;
    end
    c1 = c2;
end
w = c(:, m+1).';
end

%% =======================================================================
function [xopt, fopt] = brentMin(f, ax, bx, tol)
if nargin < 4 || isempty(tol), tol = 1e-6; end

a = min(ax,bx);
b = max(ax,bx);

CGOLD = 0.3819660112501051;
ZEPS  = 1e-12;

x = a + CGOLD*(b-a);
w = x; v = x;
fx = f(x); fw = fx; fv = fx;

d = 0.0;
e = 0.0;

for iter = 1:200
    xm = 0.5*(a+b);
    tol1 = tol*abs(x) + ZEPS;
    tol2 = 2.0*tol1;

    if abs(x - xm) <= (tol2 - 0.5*(b-a))
        break;
    end

    p = 0.0; q = 0.0; r = 0.0;
    if abs(e) > tol1
        r = (x-w)*(fx-fv);
        q = (x-v)*(fx-fw);
        p = (x-v)*q - (x-w)*r;
        q = 2.0*(q-r);
        if q > 0, p = -p; end
        q = abs(q);

        etemp = e;
        e = d;

        if (abs(p) < abs(0.5*q*etemp)) && (p > q*(a-x)) && (p < q*(b-x))
            d = p/q;
            u = x + d;
            if (u-a) < tol2 || (b-u) < tol2
                d = sign(xm-x)*tol1;
            end
        else
            if x >= xm
                e = a - x;
            else
                e = b - x;
            end
            d = CGOLD*e;
        end
    else
        if x >= xm
            e = a - x;
        else
            e = b - x;
        end
        d = CGOLD*e;
    end

    if abs(d) >= tol1
        u = x + d;
    else
        u = x + sign(d)*tol1;
    end

    fu = f(u);

    if fu <= fx
        if u >= x
            a = x;
        else
            b = x;
        end
        v = w; fv = fw;
        w = x; fw = fx;
        x = u; fx = fu;
    else
        if u < x
            a = u;
        else
            b = u;
        end
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