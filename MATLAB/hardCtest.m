function swan_optimize_refgrid_bfgs_userN
% Reference grid = [s, C t] + [sum a_k f1_k, sum b_k f2_k], with f•_k=0 on boundary.
% You choose how many basis functions: N_basis (or Kx_desired/Ky_desired).
% Outer Brent (fminbnd) optimizes C; inner custom BFGS (no toolbox) optimizes a,b.
%
% NEW: also builds & plots the REFERENCE GRID itself (Xref,Yref),
%      and prints its Jacobian stats.

clear; clc;

%% ==================== User settings ======================================
Ns = 41; Nt = 41;            % logical grid sizes (s-columns x t-rows)
N_basis = 36;                % <- total bubbles for X and for Y (set this)
Kx_desired = N_basis;        % or set these separately
Ky_desired = N_basis;

mode_order = 'radial';       % 'radial' (by m^2+n^2) or 'lexico' (by m then n)

mu = 1.0; lambda = 1.0;      % Lamé weights
C_bracket = [0.10, 2.00];    % Brent bracket for C
beta_reg = 1e-6;             % tiny Tikhonov on coeffs
epsJ = 1e-3;                 % Jacobian safety floor
kappa = 10.0;                % smooth barrier strength near J -> 0

normalize_basis = true;      % normalize bubbles by derivative RMS
print_modes = true;          % print selected & optimized modes

%% ==================== Boundaries (Swan) ==================================
I = Ns; J = Nt;
hs = 1/(I-1); ht = 1/(J-1);
nI = I-2; nJ = J-2;

xb_fun = @(s) s;                  yb_fun = @(s) 0*s;
xt_fun = @(s) s;                  yt_fun = @(s) 1 - 3*s + 3*s.^2;
xl_fun = @(s) 0*s;                yl_fun = @(s) s;
xr_fun = @(s) 1 + 2*s - 2*s.^2;   yr_fun = @(s) s;

sb = linspace(0,1,I); sl = linspace(0,1,J).';
xb = xb_fun(sb);  yb = yb_fun(sb);
xt = xt_fun(sb);  yt = yt_fun(sb);
xl = xl_fun(sl);  yl = yl_fun(sl);
xr = xr_fun(sl);  yr = yr_fun(sl);

BND.X = nan(J,I); BND.Y = nan(J,I);
BND.X(1,:) = xb;  BND.Y(1,:) = yb;
BND.X(J,:) = xt;  BND.Y(J,:) = yt;
BND.X(:,1) = xl;  BND.Y(:,1) = yl;
BND.X(:,I) = xr;  BND.Y(:,I) = yr;
BND.X(1,1)=xb(1);   BND.Y(1,1)=yb(1);
BND.X(1,I)=xb(end); BND.Y(1,I)=yb(end);
BND.X(J,1)=xt(1);   BND.Y(J,1)=yt(1);
BND.X(J,I)=xt(end); BND.Y(J,I)=yt(end);

%% ==================== Sylvester operators & RHS stencils ==================
eI = ones(nI,1); eJ = ones(nJ,1);
Ts = spdiags([-eI 2*eI -eI], [-1 0 1], nI, nI);  % s-Laplacian (Dirichlet)
Tt = spdiags([-eJ 2*eJ -eJ], [-1 0 1], nJ, nJ);  % t-Laplacian (Dirichlet)

bFxL = BND.X(2:J-1, 1);   bFyL = BND.Y(2:J-1, 1);
bFxR = BND.X(2:J-1, I);   bFyR = BND.Y(2:J-1, I);
bFxB = BND.X(1, 2:I-1);   bFyB = BND.Y(1, 2:I-1);
bFxT = BND.X(J, 2:I-1);   bFyT = BND.Y(J, 2:I-1);

%% ==================== Make bubble mode lists =============================
modes_X = make_modes(Kx_desired, mode_order);
modes_Y = make_modes(Ky_desired, mode_order);

if print_modes
    fprintf('Using %d X-basis and %d Y-basis functions.\n', size(modes_X,1), size(modes_Y,1));
    fprintf('modes_X = %s\n', mat2str(modes_X));
    fprintf('modes_Y = %s\n\n', mat2str(modes_Y));
end

%% ==================== Build bases & center-derivatives ===================
[sgrid, tgrid] = meshgrid(sb, sl);   % JxI

Kx = size(modes_X,1); Ky = size(modes_Y,1);

F1 = cell(Kx,1);  F1s_c = cell(Kx,1);  F1t_c = cell(Kx,1);
F2 = cell(Ky,1);  F2s_c = cell(Ky,1);  F2t_c = cell(Ky,1);

center_derivs = @(F) deal( ...
    (F(:,3:end)-F(:,1:end-2))/(2*hs), ...  % *_s (J x (I-2))
    (F(3:end,:)-F(1:end-2,:))/(2*ht) ...   % *_t ((J-2) x I)
);

for k = 1:Kx
    m = modes_X(k,1); n = modes_X(k,2);
    F = sin(m*pi*sgrid).*sin(n*pi*tgrid);      % zero on boundaries
    [Fs, Ft] = center_derivs(F);
    F1{k}    = F;
    F1s_c{k} = Fs(2:end-1,:);                  % (J-2) x (I-2)
    F1t_c{k} = Ft(:,2:end-1);                  % (J-2) x (I-2)
end
for k = 1:Ky
    m = modes_Y(k,1); n = modes_Y(k,2);
    F = sin(m*pi*sgrid).*sin(n*pi*tgrid);
    [Fs, Ft] = center_derivs(F);
    F2{k}    = F;
    F2s_c{k} = Fs(2:end-1,:);
    F2t_c{k} = Ft(:,2:end-1);
end

% Normalize by derivative RMS (optional)
if normalize_basis
    for k = 1:Kx
        sc = sqrt(mean(F1s_c{k}(:).^2 + F1t_c{k}(:).^2));
        if sc>0, F1{k}=F1{k}/sc; F1s_c{k}=F1s_c{k}/sc; F1t_c{k}=F1t_c{k}/sc; end
    end
    for k = 1:Ky
        sc = sqrt(mean(F2s_c{k}(:).^2 + F2t_c{k}(:).^2));
        if sc>0, F2{k}=F2{k}/sc; F2s_c{k}=F2s_c{k}/sc; F2t_c{k}=F2t_c{k}/sc; end
    end
end

cellArea = hs*ht;

%% ==================== Warm start for inner BFGS ==========================
last_z = zeros(Kx+Ky,1);   % captured by nested functions

%% ==================== Outer Brent on C; inner BFGS on (a,b) ==============
energy_at_C = @(C) inner_opt(C);   % nested (captures last_z)

opts_b = optimset('TolX',1e-3,'Display','iter');
[C_opt, E_opt] = fminbnd(@(C) energy_at_C(C), C_bracket(1), C_bracket(2), opts_b);

% Final inner solve at C* to fetch (a*,b*) and final grid
[z_star, ~] = bfgs(@(z) energy_grad_coeffs(C_opt, z), last_z, struct('display',true));
a_opt = z_star(1:Kx);  b_opt = z_star(Kx+1:end);
[~, ~, ~, Xopt, Yopt, Jc] = energy_grad_coeffs(C_opt, z_star);

fprintf('\nBest C ≈ %.6f\nTotal energy ≈ %.6e\nMin Jacobian(center) ≈ %.3e\n', C_opt, E_opt, min(Jc(:)));
fprintf('||a||_2=%.3e, ||b||_2=%.3e\n', norm(a_opt), norm(b_opt));

if print_modes
    fprintf('\nOptimized X modes [m n a]: %s\n', mat2str([modes_X, a_opt(:)], 6));
    fprintf('Optimized Y modes [m n b]: %s\n\n', mat2str([modes_Y, b_opt(:)], 6));
end

%% ==================== Plot: actual grid (physical) =======================
figure; hold on; axis equal; box on;
for j=1:J, plot(Xopt(j,:),Yopt(j,:),'-'); end
for i=1:I, plot(Xopt(:,i),Yopt(:,i),'-'); end
title(sprintf('Swan grid (physical): C*=%.3f, BFGS-optimized (N_x=%d, N_y=%d)', C_opt, Kx, Ky));
xlabel('x'); ylabel('y');

figure; imagesc(Jc); axis image; colorbar;
title('Jacobian det J at cell centers (physical grid)'); xlabel('i'); ylabel('j');

%% ==================== NEW: build & plot the REFERENCE grid ===============
[Xref, Yref] = reference_grid_from(C_opt, a_opt, b_opt);       % [s + sum a_k F1,  C t + sum b_k F2]
[~, ~, ~, ~, Jc_ref] = energy_grad_coeffs_reference(C_opt, a_opt, b_opt);

figure; hold on; axis equal; box on;
for j=1:J, plot(Xref(j,:), Yref(j,:), '-'); end
for i=1:I, plot(Xref(:,i), Yref(:,i), '-'); end
title(sprintf('REFERENCE grid: X_ref = s + Σ a_k f1_k,  Y_ref = C t + Σ b_k f2_k  (C*=%.3f)', C_opt));
xlabel('x_{ref}'); ylabel('y_{ref}');

figure; imagesc(Jc_ref); axis image; colorbar;
title('Jacobian det J at cell centers (REFERENCE grid)'); xlabel('i'); ylabel('j');

%% ==================== Nested helpers =====================================
    function modes = make_modes(K, order)
        % Generate K pairs (m,n) with m,n>=1, in chosen order.
        % 'radial': by m^2+n^2 (ties by m then n).  'lexico': by m, then n.
        maxM = ceil(2 + sqrt(K)*2);  % crude envelope
        cand = [];
        for m = 1:maxM
            for n = 1:maxM
                cand = [cand; m, n, m*m+n*n]; %#ok<AGROW>
            end
        end
        switch lower(order)
            case 'radial'
                [~,idx] = sortrows(cand,[3 1 2]); % by radius^2, then m,n
            case 'lexico'
                [~,idx] = sortrows(cand,[1 2]);   % by m, then n
            otherwise
                error('Unknown mode_order: %s', order);
        end
        modes = cand(idx(1:K),1:2);
    end

    function E = inner_opt(C)
        bfgs_opts.maxIter = 200;
        bfgs_opts.tolGrad = 1e-6;
        bfgs_opts.tolStep = 1e-10;
        bfgs_opts.display = false;
        [z_opt, E] = bfgs(@(z) energy_grad_coeffs(C, z), last_z, bfgs_opts);
        last_z = z_opt;   % warm start for next C
    end

    function [X, Y] = base_grid_from_C(C)
        wx = 1/(hs^2);             % A=1
        wy = 1/((C^2)*ht^2);
        Aop = wy * full(Tt);
        Bop = wx * full(Ts);

        Fx = zeros(nJ,nI);  Fy = zeros(nJ,nI);
        Fx(:,1)   = Fx(:,1)   + wx * bFxL;   Fy(:,1)   = Fy(:,1)   + wx * bFyL;
        Fx(:,end) = Fx(:,end) + wx * bFxR;   Fy(:,end) = Fy(:,end) + wx * bFyR;
        Fx(1,:)   = Fx(1,:)   + wy * bFxB;   Fy(1,:)   = Fy(1,:)   + wy * bFyB;
        Fx(end,:) = Fx(end,:) + wy * bFxT;   Fy(end,:) = Fy(end,:) + wy * bFyT;

        Ux = sylvester(Aop, Bop, Fx);
        Uy = sylvester(Aop, Bop, Fy);

        % Sanity: sizes must be (J-2)x(I-2)
        assert(isequal(size(Ux), [J-2, I-2]), 'Ux size mismatch');
        assert(isequal(size(Uy), [J-2, I-2]), 'Uy size mismatch');

        X = BND.X;  Y = BND.Y;
        X(2:J-1,2:I-1) = Ux;
        Y(2:J-1,2:I-1) = Uy;
    end

    function [xs, ys, xt, yt, Jc] = metrics(X, Y)
        x_s = (X(:,3:end)-X(:,1:end-2))/(2*hs);
        y_s = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
        x_t = (X(3:end,:)-X(1:end-2,:))/(2*ht);
        y_t = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);
        xs = x_s(2:end-1,:);  ys = y_s(2:end-1,:);
        xt = x_t(:,2:end-1);  yt = y_t(:,2:end-1);
        Jc = xs.*yt - xt.*ys;
    end

    % ---------- energy & gradient on the PHYSICAL grid (unchanged) ----------
    function [E, g, minJ, Xfull, Yfull, Jc] = energy_grad_coeffs(C, z)
        a = z(1:Kx); b = z(Kx+1:end);
        [X0, Y0] = base_grid_from_C(C);
        [xs, ys, xt, yt, Jc] = metrics(X0, Y0);
        for k = 1:Kx, xs = xs + a(k)*F1s_c{k}; xt = xt + a(k)*F1t_c{k}; end
        for k = 1:Ky, ys = ys + b(k)*F2s_c{k}; yt = yt + b(k)*F2t_c{k}; end
        Jsafe = max(Jc, epsJ);
        logJ  = log(Jsafe);
        d = 2;
        trC = xs.^2 + ys.^2 + xt.^2 + yt.^2;
        W = 0.5*mu*(trC - d) - mu*logJ + 0.5*lambda*(logJ.^2);
        r = max(0, epsJ - Jc);
        B = kappa * (r.^2);
        E = cellArea * sum(W(:) + B(:)) + beta_reg*(sum(a.^2)+sum(b.^2));
        minJ = min(Jc(:));
        common = (-mu + lambda*logJ) ./ Jsafe;
        dW_dxs = mu*xs + common .* (yt);
        dW_dxt = mu*xt + common .* (-ys);
        dW_dys = mu*ys + common .* (-xt);
        dW_dyt = mu*yt + common .* (xs);
        active = (Jc < epsJ);
        T = zeros(size(Jc)); T(active) = -2*kappa .* (epsJ - Jc(active));
        dW_dxs = dW_dxs + T .* (yt);
        dW_dxt = dW_dxt + T .* (-ys);
        dW_dys = dW_dys + T .* (-xt);
        dW_dyt = dW_dyt + T .* (xs);
        g = zeros(Kx+Ky,1);
        for k = 1:Kx
            g(k) = cellArea * sum( dW_dxs(:).*F1s_c{k}(:) + dW_dxt(:).*F1t_c{k}(:) ) ...
                   + 2*beta_reg*a(k);
        end
        for k = 1:Ky
            g(Kx+k) = cellArea * sum( dW_dys(:).*F2s_c{k}(:) + dW_dyt(:).*F2t_c{k}(:) ) ...
                      + 2*beta_reg*b(k);
        end
        if nargout >= 4
            Xfull = X0; Yfull = Y0;
            for k = 1:Kx, Xfull = Xfull + a(k)*F1{k}; end
            for k = 1:Ky, Yfull = Yfull + b(k)*F2{k}; end
        end
    end

    % ---------- NEW: build / evaluate the REFERENCE grid ----------
    function [Xref, Yref] = reference_grid_from(C, a, b)
        Xref = sgrid;              % s
        Yref = C * tgrid;          % C t
        for k = 1:Kx, Xref = Xref + a(k)*F1{k}; end
        for k = 1:Ky, Yref = Yref + b(k)*F2{k}; end
    end

    % convenience: compute center-Jacobian on the reference grid
    function [xs,ys,xt,yt,Jc] = metrics_reference(C,a,b)
        [Xr, Yr] = reference_grid_from(C,a,b);
        [xs,ys,xt,yt,Jc] = metrics(Xr, Yr);
    end

    % wrapper used after optimization to report J on the reference grid
    function [E, g, minJ, Xr, Yr, Jc_ref] = energy_grad_coeffs_reference(C,a,b)
        % no energy/gradient needed; just reuse code paths to get metrics
        [Xr, Yr] = reference_grid_from(C,a,b);
        [~,~,~,~,Jc_ref] = metrics(Xr,Yr);
        E = NaN; g = []; minJ = min(Jc_ref(:));
    end

    % ---------- BFGS (unchanged) ----------
    function [z, f] = bfgs(fun, z0, opts)
        if ~isfield(opts,'maxIter'),  opts.maxIter = 200; end
        if ~isfield(opts,'tolGrad'),  opts.tolGrad = 1e-6; end
        if ~isfield(opts,'tolStep'),  opts.tolStep = 1e-12; end
        if ~isfield(opts,'c1'),       opts.c1 = 1e-4; end
        if ~isfield(opts,'rho'),      opts.rho = 0.5; end
        if ~isfield(opts,'display'),  opts.display = true; end
        z = z0(:);
        [f, g] = fun(z);
        nv = numel(z);
        H  = eye(nv);
        for it = 1:opts.maxIter
            if opts.display
                fprintf('  BFGS iter %3d: f=%.6e  ||g||_inf=%.3e\n', it, f, norm(g,inf));
            end
            if norm(g,inf) < opts.tolGrad, break; end
            p = -H*g; if g'*p >= 0, p = -g; end
            t = 1.0; f0 = f; g0 = g;
            while true
                ztrial = z + t*p;
                [ftrial, gtrial] = fun(ztrial);
                if ftrial <= f0 + opts.c1*t*(g0'*p) || t < 1e-16, break; end
                t = t * opts.rho;
            end
            s = t*p;
            z = z + s;
            y = gtrial - g;
            f = ftrial; g = gtrial;
            ys = y'*s;
            if ys <= 1e-12
                H = eye(nv);
            else
                Id  = eye(nv);
                rhoB = 1/ys;
                V   = Id - rhoB*(s*y');
                H   = V*H*V' + rhoB*(s*s');
            end
            if norm(s,inf) < opts.tolStep, break; end
        end
    end
end
