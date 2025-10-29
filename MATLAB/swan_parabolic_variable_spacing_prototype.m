
function [X,Y,info] = swan_parabolic_variable_spacing_prototype(opts)
% 2‑D parabolic marching prototype (algebraic predictor) with variable spacing.
% Implements the factors εc, εi (from cumulative spacing d_k) and the
% switching εs = exp(-a*zeta). See Sec. 2.2, Eqs. (12)–(19) and (15). :contentReference[oaicite:0]{index=0}
%
% Usage:
%   [X,Y,info] = swan_parabolic_variable_spacing_prototype();
%   opts.profile = struct('type','geom','geom_r',1.08);
%   opts.a_switch = 12;  [X,Y,info] = swan_parabolic_variable_spacing_prototype(opts);

if nargin==0, opts=struct; end
I = getf(opts,'I',201);
J = getf(opts,'J',161);
a_switch = getf(opts,'a_switch',15);       % damping 'a' in εs = exp(-a*zeta)
smooth_lambda = getf(opts,'smooth_lambda',0.05);
eps0 = eps;                                % numeric safeguard

% ---- spacing profile along depth (layers 1..J; k=1 is free surface)
prof.type   = 'tanh';      % 'uniform'|'geom'|'exp'|'tanh'|'custom'
prof.alpha  = 3.0;         % tanh/exp strength
prof.geom_r = 1.06;        % geometric ratio if 'geom'
prof.custom_h = [];        % supply J-1 values if 'custom'
if isfield(opts,'profile')
    fn = fieldnames(opts.profile);
    for u = 1:numel(fn), prof.(fn{u}) = opts.profile.(fn{u}); end
end

% ---- Swan boundary curves
xb_fun = @(s) s;                  yb_fun = @(s) 0.*s;                 % bottom(s)
xt_fun = @(s) s;                  yt_fun = @(s) 1 - 3*s + 3*s.^2;     % top(s)
xl_fun = @(t) 0.*t;               yl_fun = @(t) t;                    % left(t)
xr_fun = @(t) 1 + 2*t - 2*t.^2;   yr_fun = @(t) t;                    % right(t)

% ---- logical spacings
hs = 1/(I-1); ht = 1/(J-1);
s = linspace(0,1,I);

% ---- cumulative normalized spacing d_k
[d, h] = make_profile(J, prof);          % d(1)=0, d(J)=1, sum(h)=1

% ---- build grid skeleton with d_k controlling vertical param
X = nan(J,I); Y = nan(J,I);
X(1,:) = xt_fun(s);  Y(1,:) = yt_fun(s);      % top (free surface)
X(J,:) = xb_fun(s);  Y(J,:) = yb_fun(s);      % bottom

% set sides for all layers (use t = d_k)
for k = 1:J
    t = d(k);
    X(k,1) = xl_fun(t);   Y(k,1) = yl_fun(t);
    X(k,I) = xr_fun(t);   Y(k,I) = yr_fun(t);
end

% ---- marching (algebraic predictor)
K = J;  k = 2;
while k <= K-1
    km1 = k-1;  kp1 = k+1;

    % εc, εi, εs  (Eqs. 18, 19, 15)
    eps_c = (d(k+1) - d(k-1)) / max(d(K) - d(k-1), eps0);
    eps_i = (d(k)   - d(k-1)) / max(d(k+1) - d(k-1), eps0);
    zeta  = (k-1)/(K-1);
    eps_s = exp(-a_switch * zeta);

    % data on layer k-1 and bottom K
    rkm1 = [X(km1,:).', Y(km1,:).'];     % I x 2
    rK   = [X(K,:).',   Y(K,:).'];
    RB   = rK - rkm1;                     % straight-line vector to bottom
    Rmag = sqrt(sum(RB.^2,2));            % I x 1

    % tangent along the curve at k-1 (central diff), in-plane normal
    dx = deriv_cent(X(km1,:), hs);   dy = deriv_cent(Y(km1,:), hs);
    T  = [dx(:), dy(:)];                     % I x 2
    Tn = sqrt(sum(T.^2,2));  Tn = max(Tn, 1e-14);
    N  = [-T(:,2)./Tn, T(:,1)./Tn];          % rotate +90
    % orient normal roughly toward RB
    sgn = sign(sum(N.*RB,2));  sgn(sgn==0) = 1;
    N = N .* repmat(sgn,1,2);

    % (12) orthogonal candidate: r_{k+1}^o = r_{k-1} + εc |R| n
    scale = eps_c * Rmag;                  % I x 1
    ro = rkm1 + N .* repmat(scale,1,2);    % explicit expansion (no implicit)

    % (13) straight-line candidate: r_{k+1}^s = r_{k-1} + εc R
    rs = rkm1 + eps_c * RB;

    % (14) blend with switching εs
    rkp1 = eps_s * ro + (1 - eps_s) * rs;

    % enforce side BC at layer kp1, then smooth along s
    tkp1 = d(kp1);
    rkp1(1,:)   = [xl_fun(tkp1), yl_fun(tkp1)];
    rkp1(end,:) = [xr_fun(tkp1), yr_fun(tkp1)];
    X(kp1,:) = smooth_row(rkp1(:,1).', [X(kp1,1), X(kp1,I)], smooth_lambda);
    Y(kp1,:) = smooth_row(rkp1(:,2).', [Y(kp1,1), Y(kp1,I)], smooth_lambda);

    % (16) interpolate the middle layer k between k-1 and new k+1, then smooth
    rk = rkm1 + eps_i * ([X(kp1,:).', Y(kp1,:).'] - rkm1);
    tk = d(k);
    rk(1,:)   = [xl_fun(tk), yl_fun(tk)];
    rk(end,:) = [xr_fun(tk), yr_fun(tk)];
    X(k,:) = smooth_row(rk(:,1).', [X(k,1), X(k,I)], smooth_lambda);
    Y(k,:) = smooth_row(rk(:,2).', [Y(k,1), Y(k,I)], smooth_lambda);

    k = k + 2;  % leap to the next pair
end

% if J is even, the last middle line isn't set in the loop—interpolate once
if mod(J,2) == 0
    k = J; km1 = k-1;
    eps_i = 0.5;
    rkm1 = [X(km1,:).', Y(km1,:).'];
    rk   = rkm1 + eps_i * ([X(k+1-1,:).', Y(k+1-1,:).'] - rkm1); %#ok<*NASGU>
    tk = d(k);
    rk(1,:)   = [xl_fun(tk), yl_fun(tk)];
    rk(end,:) = [xr_fun(tk), yr_fun(tk)];
    X(k,:) = smooth_row(rk(:,1).', [X(k,1), X(k,I)], smooth_lambda);
    Y(k,:) = smooth_row(rk(:,2).', [Y(k,1), Y(k,I)], smooth_lambda);
end

% ---- quality report
[Jc,~,~,~] = jacobian_from_grid(X, Y, hs, ht);
info.minJ = min(Jc(:));
info.h = h; info.d = d; info.a_switch = a_switch;

fprintf('Parabolic predictor (variable spacing): min J = %.3e\n', info.minJ);

% ---- optional plots
if nargout==0
    figure; hold on; axis equal; box on;
    for j=1:J, plot(X(j,:),Y(j,:),'-'); end
    for i=1:I, plot(X(:,i),Y(:,i),'-'); end
    title('Parabolic marching prototype with variable spacing');
    xlabel x; ylabel y;

    figure; imagesc(Jc); axis image; colorbar;
    title('Jacobian det J (cell centers)'); xlabel i; ylabel j;
end
end

% ---------------- helpers ----------------
function v = getf(s, fld, def)
if isfield(s,fld), v=s.(fld); else, v=def; end
end

function dfdx = deriv_cent(f, hs)
I = numel(f);
dfdx = zeros(1,I);
dfdx(2:I-1) = (f(3:I) - f(1:I-2))/(2*hs);
dfdx(1)     = (f(2) - f(1))/hs;
dfdx(I)     = (f(I) - f(I-1))/hs;
end

function [d,h] = make_profile(K, prof)
% returns d(1..K), h(1..K-1) with sum(h)=1, d(1)=0, d(K)=1
z = linspace(0,1,K);
switch lower(prof.type)
    case 'uniform'
        d = z;
    case 'geom'        % geometric growth from top -> bottom
        r = prof.geom_r;
        w = r.^(0:(K-2));  w = w / sum(w);
        d = [0, cumsum(w)];
    case 'exp'
        a = prof.alpha;
        d = (exp(a*z)-1)/(exp(a)-1);
    case 'tanh'
        a = prof.alpha;
        d = 0.5*(tanh(a*(2*z-1))/tanh(a) + 1);
    case 'custom'
        w = prof.custom_h(:).';
        if numel(w) ~= K-1, error('custom_h must have length K-1.'); end
        w = w / sum(w); d = [0, cumsum(w)];
    otherwise
        d = z;
end
d = d(:).'; d = d / d(end);
h = diff(d);
end

function u = smooth_row(u_pred, bnds, lam)
% implicit 1D smoothing along s: (I - lam*Dss) u = u_pred, Dirichlet ends
if lam<=0, u=u_pred; return; end
I = numel(u_pred);
u = u_pred;
u(1) = bnds(1); u(end) = bnds(2);
nI = I-2; if nI <= 0, return; end
a = -lam*ones(nI-1,1);
b = (1+2*lam)*ones(nI,1);
c = -lam*ones(nI-1,1);
d = u_pred(2:end-1).';

% forward sweep
for i=2:nI
    w = a(i-1)/b(i-1);
    b(i) = b(i) - w*c(i-1);
    d(i) = d(i) - w*d(i-1);
end
% back substitution
u_int = zeros(nI,1);
u_int(nI) = d(nI)/b(nI);
for i=nI-1:-1:1
    u_int(i) = (d(i) - c(i)*u_int(i+1))/b(i);
end
u(2:end-1) = u_int;
end

function [Jc, xs, ys, xt] = jacobian_from_grid(X, Y, hs, ht)
x_s  = (X(:,3:end)-X(:,1:end-2))/(2*hs);
y_s  = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
x_t  = (X(3:end,:)-X(1:end-2,:))/(2*ht);
y_t  = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);
xs = x_s(2:end-1,:);  ys = y_s(2:end-1,:);
xt = x_t(:,2:end-1);  yt = y_t(:,2:end-1);
Jc = xs.*yt - xt.*ys;
end