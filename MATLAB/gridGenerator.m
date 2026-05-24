function [C_opt, E_opt, Xopt, Yopt, Jc, pUopt, pVopt, Uopt, Vopt, info] = gridGenerator(grid_fun, Ns, Nt, Nsteps)
% gridGenerator
% Reference-grid reclustering driven by minimizing Var(log(J)).
%
% Main loop:
%   1) build reference grid from pU,pV
%   2) solve physical grid with anisotropic Laplace
%   3) compute Jc and Var(log(J+))
%   4) build reference-grid scores from log-J deviations
%   5) update pU,pV with backtracking
%
% The Laplace solver returns Jc and vlogJ directly.

%% -----------------------------------------------------------------------
%% initialize outputs
C_opt = NaN; E_opt = NaN;
Xopt = []; Yopt = []; Jc = [];
pUopt = []; pVopt = [];
Uopt = []; Vopt = [];
info = struct();

%% -----------------------------------------------------------------------
%% settings
I = Ns; J = Nt;
if nargin < 4 || isempty(Nsteps)
    Nsteps = 25;
end
Nsteps = max(0, floor(Nsteps));

derivOrder = 2;

thetaW = 1;
mu = 1.0; lambda = 1.0;
epsJ = 1e-6; kappa = 10.0;
epsJw = 1e-10;

gammaOrth = 4e-2;
gammaAR   = 4e-2;
betaAB    = 1e-4;

filterMode = 'minq';
qminSJ = 0.25;

epsDen = 1e-14;

hs = 1/(I-1);
ht = 1/(J-1);

useBrent = true;
C_fixed = 0.4269;
C_bracket = [0.01, 10];
C_tol = 1e-6;

% --- reference reclustering settings ---
tauJ = 1e-8;          % positivity floor for logs / rejection logic
step0 = 0.35;         % reference-step size
eqW   = 1.0;          % log-J equalization weight
badW  = 0.15;         % optional small barrier weight on bad cells
scoreSmoothPasses = 1;
maxBacktrack = 12;

doU = true;
doV = true;
% For radial/C-grid-type distortions, this is often safer:
% doU = false; doV = true;

%% -----------------------------------------------------------------------
%% boundary setup
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

%% -----------------------------------------------------------------------
%% Stage 1: choose C
pU0 = zeros(I-1,1);
pV0 = zeros(J-1,1);

if useBrent
    fC = @(C) objectiveScalar(C, pU0, pV0);
    [C_opt, ~] = brentMin(fC, C_bracket(1), C_bracket(2), C_tol);
else
    C_opt = C_fixed;
end
fprintf('C_opt = %.8f\n', C_opt);

%% -----------------------------------------------------------------------
%% Stage 2: reclustering loop
pU = pU0;
pV = pV0;

Ehist     = nan(max(Nsteps,1),1);
minJhist  = nan(max(Nsteps,1),1);
negHist   = nan(max(Nsteps,1),1);
cvJhist   = nan(max(Nsteps,1),1);
vlogJhist = nan(max(Nsteps,1),1);
stepHist  = nan(max(Nsteps,1),1);

if Nsteps <= 0
    [Uopt, Vopt] = buildReferenceGridSoftmax(C_opt, pU, pV, I, J);
    [Xopt, Yopt, Jc, ~] = solveAnisoLaplaceFromUV(Uopt, Vopt, BND, epsDen, hs, ht, derivOrder);
    pUopt = pU;
    pVopt = pV;
    E_opt = NaN;

    info.Ehist = [];
    info.minJhist = [];
    info.negHist = [];
    info.cvJhist = [];
    info.vlogJhist = [];
    info.stepHist = [];
else
    for k = 1:Nsteps
        [U, V] = buildReferenceGridSoftmax(C_opt, pU, pV, I, J);
        [X, Y, Jcur, vlogCur] = solveAnisoLaplaceFromUV(U, V, BND, epsDen, hs, ht, derivOrder);

        [Ecur, minJcur, negCur, cvCur] = evaluateStateFromJ(X, Y, Jcur, pU, pV, hs, ht, derivOrder, ...
            thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaAB);

        Ehist(k)     = Ecur;
        minJhist(k)  = minJcur;
        negHist(k)   = negCur;
        cvJhist(k)   = cvCur;
        vlogJhist(k) = vlogCur;

        [scoreUdual, scoreVdual] = logJacobianVarianceScores(Jcur, tauJ, eqW, badW);

        scoreU = dualToSpacingScore(scoreUdual);
        scoreV = dualToSpacingScore(scoreVdual);

        for s = 1:scoreSmoothPasses
            scoreU = smooth1D(scoreU);
            scoreV = smooth1D(scoreV);
        end

        if ~doU, scoreU(:) = 0; end
        if ~doV, scoreV(:) = 0; end

        qU = log(diff(U(1,:).'));
        qV = log(diff(V(:,1)));

        sc = max([max(abs(scoreU)), max(abs(scoreV)), 1e-14]);
        scoreU = scoreU / sc;
        scoreV = scoreV / sc;

        accepted = false;
        alpha = step0;

        for bt = 1:maxBacktrack
            qUtry = qU + alpha*scoreU;
            qVtry = qV + alpha*scoreV;

            qUtry = qUtry - mean(qUtry);
            qVtry = qVtry - mean(qVtry);

            [Utry, Vtry] = buildReferenceGridSoftmax(C_opt, qUtry, qVtry, I, J);
            [Xtry, Ytry, Jtry, vlogTry] = solveAnisoLaplaceFromUV(Utry, Vtry, BND, epsDen, hs, ht, derivOrder);

            [Etry, minJtry, negTry, cvTry] = evaluateStateFromJ(Xtry, Ytry, Jtry, qUtry, qVtry, hs, ht, derivOrder, ...
                thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaAB);

            if isBetterState(negTry, minJtry, vlogTry, negCur, minJcur, vlogCur)
                pU = qUtry;
                pV = qVtry;

                Ehist(k)     = Etry;
                minJhist(k)  = minJtry;
                negHist(k)   = negTry;
                cvJhist(k)   = cvTry;
                vlogJhist(k) = vlogTry;
                stepHist(k)  = alpha;

                accepted = true;
                break;
            end

            alpha = 0.5*alpha;
        end

        if ~accepted
            stepHist(k) = 0;
        end

        fprintf('iter %3d/%d: minJ=% .3e  neg=%4d  CV(J+)=% .3e  Var(logJ+)=% .3e  step=%.3e\n', ...
            k, Nsteps, minJhist(k), negHist(k), cvJhist(k), vlogJhist(k), stepHist(k));
    end

    info.Ehist     = Ehist(1:Nsteps);
    info.minJhist  = minJhist(1:Nsteps);
    info.negHist   = negHist(1:Nsteps);
    info.cvJhist   = cvJhist(1:Nsteps);
    info.vlogJhist = vlogJhist(1:Nsteps);
    info.stepHist  = stepHist(1:Nsteps);

    E_opt = Ehist(Nsteps);
    pUopt = pU;
    pVopt = pV;

    [Uopt, Vopt] = buildReferenceGridSoftmax(C_opt, pUopt, pVopt, I, J);
    [Xopt, Yopt, Jc, ~] = solveAnisoLaplaceFromUV(Uopt, Vopt, BND, epsDen, hs, ht, derivOrder);
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

%% -----------------------------------------------------------------------
%% plots
figure; hold on; axis equal; box on;
for jj=1:J, plot(Xopt(jj,:), Yopt(jj,:), '-k'); end
for ii=1:I, plot(Xopt(:,ii), Yopt(:,ii), '-k'); end
title(sprintf('Physical grid, C=%.4f', C_opt));
xlabel('X'); ylabel('Y');

figure; hold on; axis equal; box on;
for jj=1:J, plot(Uopt(jj,:), Vopt(jj,:), '-b'); end
for ii=1:I, plot(Uopt(:,ii), Vopt(:,ii), '-b'); end
title('Reference grid after log-J reclustering');
xlabel('u'); ylabel('v');

figure; imagesc(Jc); axis image; colorbar;
title('Physical Jacobian det(dX/dsdt)'); xlabel('i'); ylabel('j');

if Nsteps > 0
    figure; plot(1:Nsteps, info.minJhist, '-o'); grid on;
    xlabel('iter'); ylabel('min J');

    figure; plot(1:Nsteps, info.negHist, '-o'); grid on;
    xlabel('iter'); ylabel('# cells with J<=0');

    figure; plot(1:Nsteps, info.vlogJhist, '-o'); grid on;
    xlabel('iter'); ylabel('Var(log J_+)');

    figure; plot(1:Nsteps, info.stepHist, '-o'); grid on;
    xlabel('iter'); ylabel('accepted step');
end

%% -----------------------------------------------------------------------
%% nested objective helper
    function E = objectiveScalar(C, pUloc, pVloc)
        [U, V] = buildReferenceGridSoftmax(C, pUloc, pVloc, I, J);
        [X, Y, ~, ~] = solveAnisoLaplaceFromUV(U, V, BND, epsDen, hs, ht, derivOrder);

        ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder);
        EW  = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw);
        Ebase = (1-thetaW)*ENH + thetaW*EW;

        Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder);
        Ear   = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder);
        Ereg  = betaAB*(sum(pUloc(:).^2)+sum(pVloc(:).^2));

        E = Ebase + gammaOrth*Eorth + gammaAR*Ear + Ereg;
    end
end

%% =======================================================================
function [scoreUdual, scoreVdual] = logJacobianVarianceScores(Jc, tauJ, eqW, badW)
% Build dual-strip scores from log-J variance objective.
%
% If J is smaller than the geometric mean, score is positive:
%   allocate more reference area there.
%
% If J is larger than the geometric mean, score is negative:
%   allocate less reference area there.

Jsafe = max(Jc, tauJ);

pos = Jc > 0;
if any(pos(:))
    targetLog = mean(log(max(Jc(pos), tauJ)));
else
    targetLog = log(tauJ);
end

Req = targetLog - log(Jsafe);     % positive => cell too small
Rbad = max(0, tauJ - Jc);         % extra protection for bad cells

scoreVdual = eqW*mean(Req, 2) + badW*normalizeVec(sum(Rbad, 2));
scoreUdual = eqW*mean(Req, 1).' + badW*normalizeVec(sum(Rbad, 1).');
end

%% =======================================================================
function [E, minJ, negCount, cvPos] = evaluateStateFromJ(X, Y, Jc, pU, pV, hs, ht, derivOrder, ...
    thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaAB)

ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder);
EW  = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw);
Ebase = (1-thetaW)*ENH + thetaW*EW;

Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder);
Ear   = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder);
Ereg  = betaAB*(sum(pU(:).^2)+sum(pV(:).^2));
E = Ebase + gammaOrth*Eorth + gammaAR*Ear + Ereg;

minJ = min(Jc(:));
negCount = nnz(Jc <= 0);

pos = Jc > 0;
if any(pos(:))
    jp = Jc(pos);
    cvPos = std(jp) / max(mean(jp), 1e-14);
else
    cvPos = inf;
end
end

%% =======================================================================
function tf = isBetterState(newNeg, newMin, newVlog, oldNeg, oldMin, oldVlog)
tolMin = 1e-14;
tf = false;

if newNeg < oldNeg
    tf = true;
elseif newNeg == oldNeg && newMin > oldMin + tolMin
    tf = true;
elseif newNeg == 0 && oldNeg == 0 && newVlog < oldVlog
    tf = true;
end
end

%% =======================================================================
function s = dualToSpacingScore(d)
d = d(:);
m = numel(d) + 1;
s = zeros(m,1);

if isempty(d)
    return;
elseif numel(d) == 1
    s(:) = d(1);
    return;
end

s(1) = d(1);
s(end) = d(end);
s(2:end-1) = 0.5*(d(1:end-1) + d(2:end));
end

function v = normalizeVec(v)
v = v(:);
nrm = max(abs(v));
if nrm > 0
    v = v / nrm;
end
end

function y = smooth1D(y)
y = y(:);
if numel(y) < 3, return; end
ys = y;
ys(2:end-1) = 0.25*y(1:end-2) + 0.5*y(2:end-1) + 0.25*y(3:end);
y = ys;
end

%% =======================================================================
function [Ugrid, Vgrid] = buildReferenceGridSoftmax(C, pU, pV, I, J)
pU = pU(:);
pV = pV(:);

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
x = x - max(x);
ex = exp(x);
y = ex / sum(ex);
end

%% =======================================================================
function [minq, ESJ] = scaledJacobianStatsHO(X, Y, hs, ht, derivOrder, qmin)
if nargin < 6, qmin = 0.25; end
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
function [X, Y, Jc, vlogJ] = solveAnisoLaplaceFromUV(Ugrid, Vgrid, BND, epsDen, hs, ht, derivOrder)
% Same Laplace solve as before, but also returns:
%   Jc    = cell-centered Jacobian
%   vlogJ = variance of log(J) over positive cells

[J, I] = size(BND.X);
nI = I-2;
nJ = J-2;
N  = nI*nJ;

u = Ugrid;
v = Vgrid;

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

        ptr = ptr + 1; ii(ptr) = row; jj(ptr) = row; ss(ptr) = diagv;

        if ii_int > 1
            ptr = ptr + 1; ii(ptr) = row; jj(ptr) = idx(ii_int-1, jj_int); ss(ptr) = -Ww;
        else
            bX(row) = bX(row) + Ww * BND.X(jg,1);
            bY(row) = bY(row) + Ww * BND.Y(jg,1);
        end

        if ii_int < nI
            ptr = ptr + 1; ii(ptr) = row; jj(ptr) = idx(ii_int+1, jj_int); ss(ptr) = -We;
        else
            bX(row) = bX(row) + We * BND.X(jg,I);
            bY(row) = bY(row) + We * BND.Y(jg,I);
        end

        if jj_int > 1
            ptr = ptr + 1; ii(ptr) = row; jj(ptr) = idx(ii_int, jj_int-1); ss(ptr) = -Ws;
        else
            bX(row) = bX(row) + Ws * BND.X(1,ig);
            bY(row) = bY(row) + Ws * BND.Y(1,ig);
        end

        if jj_int < nJ
            ptr = ptr + 1; ii(ptr) = row; jj(ptr) = idx(ii_int, jj_int+1); ss(ptr) = -Wn;
        else
            bX(row) = bX(row) + Wn * BND.X(J,ig);
            bY(row) = bY(row) + Wn * BND.Y(J,ig);
        end
    end
end

Aop = sparse(ii(1:ptr), jj(1:ptr), ss(1:ptr), N, N);

Xvec = Aop \ bX;
Yvec = Aop \ bY;

X = BND.X;
Y = BND.Y;
X(2:end-1,2:end-1) = reshape(Xvec, [nI, nJ]).';
Y(2:end-1,2:end-1) = reshape(Yvec, [nI, nJ]).';

[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
Jc = xs.*yt - xt.*ys;

pos = Jc > 0;
if any(pos(:))
    vlogJ = var(log(max(Jc(pos), 1e-14)));
else
    vlogJ = inf;
end
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
r = log((a + epsr)./(b + epsr));
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
if mod(p,2) ~= 0 || p < 2
    error('derivOrder must be even >=2');
end

[J, I] = size(F);
r = p/2;
m = 2*r + 1;
Df = zeros(J, I);

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
if mod(p,2) ~= 0 || p < 2
    error('derivOrder must be even >=2');
end

[J, I] = size(F);
r = p/2;
m = 2*r + 1;
Df = zeros(J, I);

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
if nargin < 4 || isempty(tol)
    tol = 1e-6;
end

a = min(ax,bx);
b = max(ax,bx);

CGOLD = 0.3819660112501051;
ZEPS  = 1e-12;

x = a + CGOLD*(b-a);
w = x;
v = x;
fx = f(x);
fw = fx;
fv = fx;

d = 0.0;
e = 0.0;

for iter = 1:200
    xm = 0.5*(a+b);
    tol1 = tol*abs(x) + ZEPS;
    tol2 = 2.0*tol1;

    if abs(x - xm) <= (tol2 - 0.5*(b-a))
        break;
    end

    p = 0.0;
    q = 0.0;
    r = 0.0;

    if abs(e) > tol1
        r = (x-w)*(fx-fv);
        q = (x-v)*(fx-fw);
        p = (x-v)*q - (x-w)*r;
        q = 2.0*(q-r);
        if q > 0
            p = -p;
        end
        q = abs(q);

        etemp = e;
        e = d;

        if (abs(p) < abs(0.5*q*etemp)) && (p > q*(a-x)) && (p < q*(b-x))
            d = p/q;
            u = x + d;
            if (u-a) < tol2 || (b-u) < tol2
                d = sign_nonzero(xm-x)*tol1;
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
        u = x + sign_nonzero(d)*tol1;
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

function s = sign_nonzero(x)
if x >= 0
    s = 1;
else
    s = -1;
end
end