function swan_optimize_B_bubbles_trisection
% Bubble-augmented reference grid optimization with feasibility-aware TRISECTION over B.
% Base reference = [s, B t]. Bubbles: X += sum_k a_k f_k, Y += sum_k b_k f_k,
% f_k = sin(m*pi*s) sin(n*pi*t) (zero on the boundary).
% Inner: toolbox-free BFGS over z=[a;b]. Negative Jacobians => large finite penalty (no Inf/NaN).
% Outer: trisection looks for a feasible window then shrinks using ternary logic, avoiding penalized sides.

clear; clc;

%% ---------------- User settings ----------------
I = 11; J = 11;                    % logical grid size (s-columns x t-rows)
A = 1.0;                           % fixed, only B varies
B_bracket = [0.10, 2.00];          % search range for B

mu = 1.0; lambda = 1.0;            % Lame weights for neo-Hookean energy

N_basis = 0;                      % number of bubble modes for X and for Y (same set)
normalize_basis = true;            % scale bubbles by derivative RMS for stability
print_modes = true;

% Penalty for any negative Jacobian (centers or discrete cells)
PENALTY_BAD = 1e6;                 % base penalty
NEG_COUNT_WEIGHT = 1e3;            % extra per bad location

% Small-positive cell barrier (encourages margin from zero area)
use_orient_barrier = true; eps_orient = 1e-3; kpos = 1e2; ppos = 4;

% Regularization on bubble coefficients
beta_reg = 1e-8;

% Trisection (ternary search) settings
tri.tolB   = 1e-3;
tri.maxIt  = 60;
tri.Nscan  = 33;

%% ---------------- Boundaries (Swan) ----------------
hs = 1/(I-1); ht = 1/(J-1);
nI = I-2; nJ = J-2;

xb_fun = @(s) s;                  yb_fun = @(s) 0*s;                 % bottom
xt_fun = @(s) s;                  yt_fun = @(s) 1 - 3*s + 3*s.^2;    % top
xl_fun = @(t) 0*t;                yl_fun = @(t) t;                   % left
xr_fun = @(t) 1 + 2*t - 2*t.^2;   yr_fun = @(t) t;                   % right

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

%% ---------------- Sylvester operators ----------------
eI = ones(nI,1); eJ = ones(nJ,1);
Ts = spdiags([-eI 2*eI -eI], [-1 0 1], nI, nI);
Tt = spdiags([-eJ 2*eJ -eJ], [-1 0 1], nJ, nJ);

bFxL = BND.X(2:J-1, 1);   bFyL = BND.Y(2:J-1, 1);
bFxR = BND.X(2:J-1, I);   bFyR = BND.Y(2:J-1, I);
bFxB = BND.X(1, 2:I-1);   bFyB = BND.Y(1, 2:I-1);
bFxT = BND.X(J, 2:I-1);   bFyT = BND.Y(J, 2:I-1);

%% ---------------- Bubble basis (zero on boundary) ----------------
[sgrid, tgrid] = meshgrid(sb, sl);
modes = make_modes(N_basis, 'radial');   % pairs (m,n), same set for X and Y
K = size(modes,1);

F = cell(K,1);  Fs_c = cell(K,1);  Ft_c = cell(K,1);
for k = 1:K
    m = modes(k,1); n = modes(k,2);
    Fi = sin(m*pi*sgrid).*sin(n*pi*tgrid);         % vanishes on boundary
    [Fs, Ft] = center_derivs(Fi, hs, ht);          % Jx(I-2) and (J-2)xI
    F{k}    = Fi;
    Fs_c{k} = Fs(2:end-1,:);                       % (J-2)x(I-2) at cell centers
    Ft_c{k} = Ft(:,2:end-1);                       % (J-2)x(I-2) at cell centers
end

if normalize_basis
    for k = 1:K
        sc = sqrt(mean(Fs_c{k}(:).^2 + Ft_c{k}(:).^2));
        if sc>0
            F{k} = F{k}/sc; Fs_c{k} = Fs_c{k}/sc; Ft_c{k} = Ft_c{k}/sc;
        end
    end
end
if print_modes
    fprintf('Using %d bubbles for X and %d for Y\n', K, K);
    fprintf('modes = %s\n\n', mat2str(modes));
end

%% ---------------- Pack context for local functions ----------------
ctx.I = I; ctx.J = J; ctx.hs = hs; ctx.ht = ht;
ctx.nI = nI; ctx.nJ = nJ;
ctx.BND = BND;
ctx.Ts = Ts; ctx.Tt = Tt;
ctx.bFxL = bFxL; ctx.bFyL = bFyL; ctx.bFxR = bFxR; ctx.bFyR = bFyR;
ctx.bFxB = bFxB; ctx.bFyB = bFyB; ctx.bFxT = bFxT; ctx.bFyT = bFyT;
ctx.mu = mu; ctx.lambda = lambda;
ctx.A = A;
ctx.K = K;
ctx.F = F; ctx.Fs_c = Fs_c; ctx.Ft_c = Ft_c;
ctx.PENALTY_BAD = PENALTY_BAD;
ctx.NEG_COUNT_WEIGHT = NEG_COUNT_WEIGHT;
ctx.use_orient_barrier = use_orient_barrier;
ctx.eps_orient = eps_orient; ctx.kpos = kpos; ctx.ppos = ppos;
ctx.beta_reg = beta_reg;

% Initialize persistent storage for inner solver
inner_min_for_B('init', ctx);

%% ---------------- Trisection over B with feasibility awareness -------------
[B_opt, E_opt, info] = trisection_min_feasible(@(B) funE_returnE(@inner_min_for_B, B), ...
                                               B_bracket(1), B_bracket(2), tri, PENALTY_BAD);

if ~info.found_feasible
    warning('No feasible (non-penalized) B found in [%.3f, %.3f]. Returning best scanned point.', ...
            B_bracket(1), B_bracket(2));
end
% --- Add rectangular reference grid [s, B*t] overlay ----------------------
[Sg, Tg] = meshgrid(sb, sl);     % sb: 1xI, sl: Jx1 (already defined above)
Xref = Sg; 
Yref = B_opt * Tg;

figure; hold on; axis tight; box on;
% reference grid (dashed, light)
for j0 = 1:J, plot(Xref(j0,:), Yref(j0,:), '-', 'LineWidth', 0.6); end
for i0 = 1:I, plot(Xref(:,i0), Yref(:,i0), '-', 'LineWidth', 0.6); end
title(sprintf('Best rectangular reference grid')); xlabel('x'); ylabel('y');


% Build final grid for best B and print metrics
[~, ~, Xopt, Yopt, Jc] = inner_min_for_B('eval_full', B_opt);
fprintf('\nBest B ~= %.6f\nTotal energy ~= %.6e\nMin Jacobian at centers ~= %.3e\n', ...
        B_opt, E_opt, min(Jc(:)));
% Plots
figure; hold on; axis equal; box on;
for j0=1:J, plot(Xopt(j0,:),Yopt(j0,:),'-'); end
for i0=1:I, plot(Xopt(:,i0),Yopt(:,i0),'-'); end
title(sprintf('Swan grid with bubbles at B*=%.3f (A=1).', B_opt)); xlabel('x'); ylabel('y');

figure; imagesc(Jc); axis image; colorbar;
title('Jacobian det J at cell centers'); xlabel('i'); ylabel('j');

end  % ===== end main =======================================================


%% ======================= Local functions =================================

function modesMN = make_modes(Kin, order)
    maxM = max(3, ceil(2 + sqrt(Kin)*2));
    cand = [];
    for m=1:maxM
        for n=1:maxM
            cand = [cand; m, n, m*m+n*n]; %#ok<AGROW>
        end
    end
    switch lower(order)
        case 'radial', [~,idx] = sortrows(cand,[3 1 2]);
        case 'lexico', [~,idx] = sortrows(cand,[1 2]);
        otherwise, error('Unknown order');
    end
    modesMN = cand(idx(1:Kin),1:2);
end

function [Fs, Ft] = center_derivs(U, hs_, ht_)
    Fs = (U(:,3:end)-U(:,1:end-2))/(2*hs_);   % J x (I-2)
    Ft = (U(3:end,:)-U(1:end-2,:))/(2*ht_);   % (J-2) x I
end

function [xs, ys, xt, yt, Jc] = metrics_center(X, Y, hs, ht)
    x_s = (X(:,3:end)-X(:,1:end-2))/(2*hs);
    y_s = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
    x_t = (X(3:end,:)-X(1:end-2,:))/(2*ht);
    y_t = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);
    xs = x_s(2:end-1,:);  ys = y_s(2:end-1,:);
    xt = x_t(:,2:end-1);  yt = y_t(:,2:end-1);
    Jc = xs.*yt - xt.*ys;
end

function [Jcell_norm, Ex_ll, Ey_ll, Fx_ll, Fy_ll] = cell_orientation_LL(X, Y, hs, ht)
    [J, I] = size(X);
    Ex = X(:,2:end) - X(:,1:end-1);   % J x (I-1)
    Ey = Y(:,2:end) - Y(:,1:end-1);
    Fx = X(2:end,:) - X(1:end-1,:);   % (J-1) x I
    Fy = Y(2:end,:) - Y(1:end-1,:);
    Ex_ll = Ex(1:J-1, 1:I-1); Ey_ll = Ey(1:J-1, 1:I-1);
    Fx_ll = Fx(1:J-1, 1:I-1); Fy_ll = Fy(1:J-1, 1:I-1);
    Jcell_area = Ex_ll .* Fy_ll - Ey_ll .* Fx_ll;
    Jcell_norm = Jcell_area / (hs*ht);
end

function outE = funE_returnE(funInnerB, B)
    % funInnerB is the inner_min_for_B dispatcher below
    [E, ~, ~, ~, ~] = funInnerB('eval', B);
    outE = E;
end

function [B_best, E_best, info] = trisection_min_feasible(funE, a, b, tri, PENALTY_BAD)
    info = struct('found_feasible', false, 'iterations', 0);

    % scan
    Bs = linspace(a, b, tri.Nscan);
    Es = zeros(size(Bs));
    mask = false(size(Bs));
    for i = 1:numel(Bs)
        Ei = funE(Bs(i));
        Es(i) = Ei;
        mask(i) = (Ei < PENALTY_BAD);
    end

    if any(mask)
        blocks = contiguous_blocks(mask);
        best_block = [];
        best_minE = inf;
        for k = 1:size(blocks,1)
            idx1 = blocks(k,1); idx2 = blocks(k,2);
            localMin = min(Es(idx1:idx2));
            if localMin < best_minE
                best_minE = localMin; best_block = [idx1 idx2];
            end
        end
        A = Bs(best_block(1)); B = Bs(best_block(2));
        info.found_feasible = true;

        for it = 1:tri.maxIt
            info.iterations = it;
            if (B - A) <= tri.tolB, break; end
            m1 = A + (B-A)/3;
            m2 = B - (B-A)/3;

            E1 = funE(m1);
            E2 = funE(m2);

            f1_bad = (E1 >= PENALTY_BAD);
            f2_bad = (E2 >= PENALTY_BAD);

            if ~f1_bad && ~f2_bad
                if E1 < E2
                    B = m2;
                else
                    A = m1;
                end
            elseif f1_bad && ~f2_bad
                A = m1;      % push right, left interior is bad
            elseif ~f1_bad && f2_bad
                B = m2;      % push left, right interior is bad
            else
                mc = 0.5*(A+B);
                Ec = funE(mc);
                if Ec < PENALTY_BAD
                    A = A + 0.25*(mc - A);
                    B = B - 0.25*(B - mc);
                else
                    mid = 0.5*(A+B);
                    A = A + 0.25*(B - A);
                    B = mid + 0.25*(B - mid);
                end
            end
        end

        B_best = 0.5*(A+B);
        E_best = funE(B_best);
        if E_best >= PENALTY_BAD
            feasIdx = find(mask);
            [Emin, relIdx] = min(Es(mask));
            B_best = Bs(feasIdx(relIdx));
            E_best = Emin;
        end
    else
        [E_best, idxMin] = min(Es);
        B_best = Bs(idxMin);
    end
end

function blocks = contiguous_blocks(mask)
    mask = mask(:).';
    d = diff([false, mask, false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    blocks = [starts(:), ends(:)];
end

% ====================== Inner solver dispatcher ===========================
function varargout = inner_min_for_B(mode, arg)
    % Persistent store for context and warm-start
    persistent ctx last_z

    switch mode
        case 'init'
            ctx = arg;                        % struct with all data
            last_z = zeros(2*ctx.K,1);
            varargout = {[]};

        case 'eval'
            B = arg;
            [E, minJ, X, Y, Jc, z_out] = inner_min_eval(ctx, last_z, B);
            last_z = z_out;                  % warm start for next B
            varargout = {E, minJ, X, Y, Jc};

        case 'eval_full'
            B = arg;
            [E, minJ, X, Y, Jc, ~] = inner_min_eval(ctx, last_z, B);
            varargout = {E, minJ, X, Y, Jc};

        otherwise
            error('inner_min_for_B: unknown mode');
    end
end


function [E, minJ, X, Y, Jc, z_opt] = inner_min_eval(ctx, last_z, B)
    % Build base grid by Sylvester for this B
    wx = 1/(ctx.A^2*ctx.hs^2);
    wy = 1/((B^2)*ctx.ht^2);
    Aop = wy * full(ctx.Tt);
    Bop = wx * full(ctx.Ts);

    X = ctx.BND.X; Y = ctx.BND.Y;
    Fx = zeros(ctx.nJ,ctx.nI);  Fy = zeros(ctx.nJ,ctx.nI);
    Fx(:,1)   = Fx(:,1)   + wx * ctx.bFxL;   Fy(:,1)   = Fy(:,1)   + wx * ctx.bFyL;
    Fx(:,end) = Fx(:,end) + wx * ctx.bFxR;   Fy(:,end) = Fy(:,end) + wx * ctx.bFyR;
    Fx(1,:)   = Fx(1,:)   + wy * ctx.bFxB;   Fy(1,:)   = Fy(1,:)   + wy * ctx.bFyB;
    Fx(end,:) = Fx(end,:) + wy * ctx.bFxT;   Fy(end,:) = Fy(end,:) + wy * ctx.bFyT;
    Ux = sylvester(Aop, Bop, Fx);
    Uy = sylvester(Aop, Bop, Fy);
    X(2:ctx.J-1,2:ctx.I-1) = Ux;
    Y(2:ctx.J-1,2:ctx.I-1) = Uy;

    % Optimize bubbles (toolbox-free BFGS)
    opts.maxIter=200; opts.tolGrad=1e-6; opts.tolStep=1e-12; opts.display=false;
    base.X0 = X; base.Y0 = Y;
    [z_opt, E] = bfgs(@(z) energy_grad(ctx, base, z), last_z, opts);

    % Return final grid and Jacobian at centers
    [~, ~, X, Y, Jc] = energy_grad(ctx, base, z_opt);
    minJ = min(Jc(:));
end

function [E, g, Xfull, Yfull, Jc] = energy_grad(ctx, base, z)
    % z = [a(1:K); b(1:K)]
    K = ctx.K;
    a = z(1:K); b = z(K+1:end);

    % Center metrics (base)
    [xs, ys, xt, yt, ~] = metrics_center(base.X0, base.Y0, ctx.hs, ctx.ht);

    % Add bubble derivatives at centers
    xs_ = xs; xt_ = xt; ys_ = ys; yt_ = yt;
    for kk=1:K, xs_ = xs_ + a(kk)*ctx.Fs_c{kk}; xt_ = xt_ + a(kk)*ctx.Ft_c{kk}; end
    for kk=1:K, ys_ = ys_ + b(kk)*ctx.Fs_c{kk}; yt_ = yt_ + b(kk)*ctx.Ft_c{kk}; end

    % Jacobian at centers
    Jc = xs_.*yt_ - xt_.*ys_;

    % Build full nodal fields for discrete cell orientation (LL indexing)
    Xfull = base.X0; Yfull = base.Y0;
    for kk=1:K, Xfull = Xfull + a(kk)*ctx.F{kk}; end
    for kk=1:K, Yfull = Yfull + b(kk)*ctx.F{kk}; end

    [Jcell_norm, Ex_ll, Ey_ll, Fx_ll, Fy_ll] = cell_orientation_LL(Xfull, Yfull, ctx.hs, ctx.ht);
    neg_cnt = nnz(Jc(:) <= 0) + nnz(Jcell_norm(:) <= 0);

    % Penalize if any negative Jacobian
    if neg_cnt > 0
        E = ctx.PENALTY_BAD + ctx.NEG_COUNT_WEIGHT * double(neg_cnt);
        g = zeros(2*K,1);
        return;
    end

    % Neo-Hookean energy at centers (safe: all J>0)
    trC = xs_.^2 + ys_.^2 + xt_.^2 + yt_.^2;
    logJ = log(Jc);
    d=2;
    W = 0.5*ctx.mu*(trC - d) - ctx.mu*logJ + 0.5*ctx.lambda*(logJ.^2);
    cellArea = ctx.hs*ctx.ht;
    E_soft = cellArea * sum(W(:));

    % Gradient wrt metrics
    common = (-ctx.mu + ctx.lambda*logJ) ./ Jc;
    dW_dxs = ctx.mu*xs_ + common .* (yt_);
    dW_dxt = ctx.mu*xt_ + common .* (-ys_);
    dW_dys = ctx.mu*ys_ + common .* (-xt_);
    dW_dyt = ctx.mu*yt_ + common .* (xs_);

    % Assemble gradient wrt a,b from the center derivatives of bubbles
    g = zeros(2*K,1);
    for kk=1:K
        g(kk)     = cellArea * sum( dW_dxs(:).*ctx.Fs_c{kk}(:) + dW_dxt(:).*ctx.Ft_c{kk}(:) );
        g(K+kk)   = cellArea * sum( dW_dys(:).*ctx.Fs_c{kk}(:) + dW_dyt(:).*ctx.Ft_c{kk}(:) );
    end

    % Small-positive orientation barrier (push away from ~0)
    if ctx.use_orient_barrier
        mask_pos = (Jcell_norm > 0) & (Jcell_norm < ctx.eps_orient);
        if any(mask_pos(:))
            dpos = (ctx.eps_orient - Jcell_norm(mask_pos));
            E_bar = sum( ctx.kpos * (dpos.^ctx.ppos) );

            % Chain rule: dJcell/da, dJcell/db using LL edges
            dB_dJ = zeros(size(Jcell_norm));
            dB_dJ(mask_pos) = -ctx.kpos * ctx.ppos * (dpos.^(ctx.ppos-1));
            scaleA = 1/(ctx.hs*ctx.ht);

            for kk=1:K
                dEx_ll = (ctx.F{kk}(:,2:end) - ctx.F{kk}(:,1:end-1)); dEx_ll = dEx_ll(1:ctx.J-1, 1:ctx.I-1);
                dFx_ll = (ctx.F{kk}(2:end,:) - ctx.F{kk}(1:end-1,:)); dFx_ll = dFx_ll(1:ctx.J-1, 1:ctx.I-1);
                dJ_da = dEx_ll .* Fy_ll - Ey_ll .* dFx_ll;
                g(kk)  = g(kk) + scaleA * sum( dB_dJ(:) .* dJ_da(:) );
            end
            for kk=1:K
                dEy_ll = (ctx.F{kk}(:,2:end) - ctx.F{kk}(:,1:end-1)); dEy_ll = dEy_ll(1:ctx.J-1, 1:ctx.I-1);
                dFy_ll = (ctx.F{kk}(2:end,:) - ctx.F{kk}(1:end-1,:)); dFy_ll = dFy_ll(1:ctx.J-1, 1:ctx.I-1);
                dJ_db = Ex_ll .* dFy_ll - dEy_ll .* Fx_ll;
                g(K+kk) = g(K+kk) + scaleA * sum( dB_dJ(:) .* dJ_db(:) );
            end

            E_soft = E_soft + E_bar;
        end
    end

    % Regularization
    E_reg = ctx.beta_reg*(sum(a.^2)+sum(b.^2));
    g(1:K)       = g(1:K)       + 2*ctx.beta_reg*a;
    g(K+1:end)   = g(K+1:end)   + 2*ctx.beta_reg*b;

    E = E_soft + E_reg;
end

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
            fprintf('  BFGS %3d: f=%.6e  ||g||_inf=%.3e\n', it, f, norm(g,inf));
        end
        if ~isfinite(f), break; end
        if norm(g,inf) < opts.tolGrad, break; end

        p = -H*g; if g'*p >= 0, p = -g; end

        % Armijo backtracking (finite penalties handled naturally)
        t = 1.0; f0 = f; g0 = g;
        while true
            ztrial = z + t*p;
            [ftrial, gtrial] = fun(ztrial);
            if isfinite(ftrial) && (ftrial <= f0 + opts.c1*t*(g0'*p)), break; end
            t = t * opts.rho;
            if t < 1e-16, break; end
        end

        s = t*p; z = z + s;
        y = gtrial - g;
        f = ftrial; g = gtrial;

        ys = y'*s;
        if ys <= 1e-12
            H = eye(nv);
        else
            Id = eye(nv);
            rhoB = 1/ys;
            V = Id - rhoB*(s*y');
            H = V*H*V' + rhoB*(s*s');
        end

        if norm(s,inf) < opts.tolStep, break; end
    end
end

