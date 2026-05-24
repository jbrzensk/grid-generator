function [A_opt, C_opt, E_opt, Xopt, Yopt, Jc, pUopt, pVopt, pBopt, pTopt, pLopt, pRopt, Uopt, Vopt, info] = gridGeneratorMultigrid(grid_fun, Ns, Nt, Nsteps)
% gridGeneratorMultigrid
% Simultaneous reference-grid and physical-boundary-grid Jacobian equalization.
%
% Compared with the original gridGenerator, this version has two coupled grids:
%
%   1) Reference grid:
%        U,V are generated from log-spacing variables pU,pV.
%
%   2) Physical boundary grid:
%        the four physical boundary curves are reparameterized by log-spacing
%        variables pB,pT,pL,pR. The PDE then fills the interior.
%
% The SGD/backtracking equalization step updates all six spacing variables at
% once:
%
%        pU, pV, pB, pT, pL, pR.
%
% The physical interior is still obtained from the anisotropic Dirichlet solve
%
%        X_ss + (1/A^2) X_tt = 0,
%        Y_ss + (1/C^2) Y_tt = 0.
%
% The script first chooses A,C by minimizing the regularized Neo-Hookean
% objective on the initial grids, then holds A,C fixed during multigrid
% Jacobian equalization.

%% -----------------------------------------------------------------------
%% initialize outputs
A_opt = NaN; C_opt = NaN; E_opt = NaN;
Xopt = []; Yopt = []; Jc = [];
pUopt = []; pVopt = [];
pBopt = []; pTopt = []; pLopt = []; pRopt = [];
Uopt = []; Vopt = [];
info = struct();

%% -----------------------------------------------------------------------
%% settings 
I = Ns;
J = Nt;

if nargin < 4 || isempty(Nsteps)
    Nsteps = 25;
end
Nsteps = max(0, floor(Nsteps));

derivOrder = 2;

% Energy settings
mu = 1.0;
lambda = 1.0;
epsJ = 1e-6;
kappa = 10.0;
epsJw = 1e-10;
thetaW = 0.5;              % 0 = Neo-Hookean, 1 = Winslow

gammaOrth = 4e-2;
gammaAR   = 4e-2;

% Regularization for the six log-spacing vectors.
betaRef  = 1e-4;
betaPhys = 1e-4;

filterMode = 'minq';
qminSJ = 0.25;

epsDen = 1e-14;
hs = 1/(I-1);
ht = 1/(J-1);

% PDE anisotropy optimization.
useACSearch = true;
A_fixed = 1.5006;
C_fixed = 0.1035;
A_bracket = [0.05, 20];
C_bracket = [0.05, 20];
nAcoarse = 20;
nCcoarse = 20;
AC_fmin_tol = 1e-5;

% Reference-grid height. Keep this separate from the PDE coefficient C.
refHeight = 1.0;

% Multigrid equalization settings.
tauJ = 1e-8;
step0 = 0.35;
eqW = 1.0;
badW = 0.15;
scoreSmoothPasses = 1;
maxBacktrack = 12;

% Turn individual updates on/off.
doRefU = true;
doRefV = true;
doPhysBottom = true;
doPhysTop    = true;
doPhysLeft   = true;
doPhysRight  = true;

%% -----------------------------------------------------------------------
%% boundary setup
[Xc, Yc] = grid_fun(I-1, J-1);

if ~isequal(size(Xc), [J I]) && isequal(size(Xc), [I J])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.';
    Yc = Yc.';
end

if ~isequal(size(Xc), [J I])
    error('grid_fun must return size (Nt x Ns) = (%d x %d). Got %s.', J, I, mat2str(size(Xc)));
end

% Store the original four physical boundary curves. These curves are fixed as
% geometric curves, but their nodal distributions are optimized.
CURVES.bottomX = Xc(1,:).';      CURVES.bottomY = Yc(1,:).';
CURVES.topX    = Xc(end,:).';    CURVES.topY    = Yc(end,:).';
CURVES.leftX   = Xc(:,1);        CURVES.leftY   = Yc(:,1);
CURVES.rightX  = Xc(:,end);      CURVES.rightY  = Yc(:,end);

% Exact corner values.
CURVES.cBL = [Xc(1,1),     Yc(1,1)];
CURVES.cBR = [Xc(1,end),   Yc(1,end)];
CURVES.cTL = [Xc(end,1),   Yc(end,1)];
CURVES.cTR = [Xc(end,end), Yc(end,end)];

%% -----------------------------------------------------------------------
%% Initial spacing variables
pU0 = zeros(I-1,1);
pV0 = zeros(J-1,1);
pB0 = zeros(I-1,1);
pT0 = zeros(I-1,1);
pL0 = zeros(J-1,1);
pR0 = zeros(J-1,1);

%% -----------------------------------------------------------------------
%% Stage 1: choose A,C on the initial coupled grids
if useACSearch
    Avec = logspace(log10(A_bracket(1)), log10(A_bracket(2)), nAcoarse);
    Cvec = logspace(log10(C_bracket(1)), log10(C_bracket(2)), nCcoarse);

    bestE = inf;
    bestA = A_fixed;
    bestC = C_fixed;

    fprintf('Coarse search over A,C...\n');
    for ia = 1:numel(Avec)
        for ic = 1:numel(Cvec)
            Atest = Avec(ia);
            Ctest = Cvec(ic);
            Etest = objectiveAC(Atest, Ctest, pU0, pV0, pB0, pT0, pL0, pR0);
            if Etest < bestE
                bestE = Etest;
                bestA = Atest;
                bestC = Ctest;
            end
        end
    end

    fprintf('Coarse best: A = %.8f, C = %.8f, E = %.8e\n', bestA, bestC, bestE);

    z0 = log([bestA; bestC]);
    opts = optimset('Display', 'iter', ...
                    'TolX', AC_fmin_tol, ...
                    'TolFun', AC_fmin_tol, ...
                    'MaxIter', 100, ...
                    'MaxFunEvals', 300);
    zopt = fminsearch(@(z) objectiveLogAC(z, pU0, pV0, pB0, pT0, pL0, pR0), z0, opts);

    A_opt = exp(zopt(1));
    C_opt = exp(zopt(2));
else
    A_opt = A_fixed;
    C_opt = C_fixed;
end

fprintf('A_opt = %.8f, C_opt = %.8f\n', A_opt, C_opt);

%% -----------------------------------------------------------------------
%% Stage 2: simultaneous reference-grid and physical-grid equalization
pU = pU0;
pV = pV0;
pB = pB0;
pT = pT0;
pL = pL0;
pR = pR0;

Ehist     = nan(max(Nsteps,1),1);
minJhist  = nan(max(Nsteps,1),1);
negHist   = nan(max(Nsteps,1),1);
cvJhist   = nan(max(Nsteps,1),1);
vlogJhist = nan(max(Nsteps,1),1);
stepHist  = nan(max(Nsteps,1),1);

if Nsteps <= 0
    [Uopt, Vopt] = buildReferenceGridSoftmax(refHeight, pU, pV, I, J);
    BND = buildPhysicalBoundaryFromSoftmax(CURVES, pB, pT, pL, pR, I, J);
    [Xopt, Yopt, Jc, ~] = solveAnisoLaplaceACFromUV(Uopt, Vopt, BND, A_opt, C_opt, epsDen, hs, ht, derivOrder);
    E_opt = objectiveAC(A_opt, C_opt, pU, pV, pB, pT, pL, pR);
else
    for k = 1:Nsteps
        [U, V] = buildReferenceGridSoftmax(refHeight, pU, pV, I, J);
        BND = buildPhysicalBoundaryFromSoftmax(CURVES, pB, pT, pL, pR, I, J);
        [X, Y, Jcur, vlogCur] = solveAnisoLaplaceACFromUV(U, V, BND, A_opt, C_opt, epsDen, hs, ht, derivOrder);

        [Ecur, minJcur, negCur, cvCur] = evaluateStateFromJ(X, Y, Jcur, pU, pV, pB, pT, pL, pR, hs, ht, derivOrder, ...
            thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaRef, betaPhys);

        Ehist(k)     = Ecur;
        minJhist(k)  = minJcur;
        negHist(k)   = negCur;
        cvJhist(k)   = cvCur;
        vlogJhist(k) = vlogCur;

        % Cell-centered log-J residuals generate dual-strip scores. These are
        % then converted to nodal spacing scores.
        [scoreUdual, scoreVdual] = logJacobianVarianceScores(Jcur, tauJ, eqW, badW);

        scoreRefU = dualToSpacingScore(scoreUdual);
        scoreRefV = dualToSpacingScore(scoreVdual);

        % Multigrid part:
        % Use the same column/row pressure to reparameterize the physical
        % boundary curves. Bottom/top receive horizontal scores; left/right
        % receive vertical scores.
        scoreB = scoreRefU;
        scoreT = scoreRefU;
        scoreL = scoreRefV;
        scoreR = scoreRefV;

        for s = 1:scoreSmoothPasses
            scoreRefU = smooth1D(scoreRefU);
            scoreRefV = smooth1D(scoreRefV);
            scoreB = smooth1D(scoreB);
            scoreT = smooth1D(scoreT);
            scoreL = smooth1D(scoreL);
            scoreR = smooth1D(scoreR);
        end

        if ~doRefU,       scoreRefU(:) = 0; end
        if ~doRefV,       scoreRefV(:) = 0; end
        if ~doPhysBottom, scoreB(:)    = 0; end
        if ~doPhysTop,    scoreT(:)    = 0; end
        if ~doPhysLeft,   scoreL(:)    = 0; end
        if ~doPhysRight,  scoreR(:)    = 0; end

        sc = max([max(abs(scoreRefU)), max(abs(scoreRefV)), ...
                  max(abs(scoreB)),    max(abs(scoreT)), ...
                  max(abs(scoreL)),    max(abs(scoreR)), 1e-14]);

        scoreRefU = scoreRefU / sc;
        scoreRefV = scoreRefV / sc;
        scoreB    = scoreB    / sc;
        scoreT    = scoreT    / sc;
        scoreL    = scoreL    / sc;
        scoreR    = scoreR    / sc;

        accepted = false;
        alpha = step0;

        for bt = 1:maxBacktrack
            qU = normalizeLogSpacing(pU + alpha*scoreRefU);
            qV = normalizeLogSpacing(pV + alpha*scoreRefV);
            qB = normalizeLogSpacing(pB + alpha*scoreB);
            qT = normalizeLogSpacing(pT + alpha*scoreT);
            qL = normalizeLogSpacing(pL + alpha*scoreL);
            qR = normalizeLogSpacing(pR + alpha*scoreR);

            [Utry, Vtry] = buildReferenceGridSoftmax(refHeight, qU, qV, I, J);
            BNDtry = buildPhysicalBoundaryFromSoftmax(CURVES, qB, qT, qL, qR, I, J);
            [Xtry, Ytry, Jtry, vlogTry] = solveAnisoLaplaceACFromUV(Utry, Vtry, BNDtry, A_opt, C_opt, epsDen, hs, ht, derivOrder);

            [Etry, minJtry, negTry, cvTry] = evaluateStateFromJ(Xtry, Ytry, Jtry, qU, qV, qB, qT, qL, qR, hs, ht, derivOrder, ...
                thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaRef, betaPhys);

            if isBetterState(negTry, minJtry, vlogTry, Etry, negCur, minJcur, vlogCur, Ecur)
                pU = qU;
                pV = qV;
                pB = qB;
                pT = qT;
                pL = qL;
                pR = qR;

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

        fprintf('iter %3d/%d: E=% .6e  minJ=% .3e  neg=%4d  CV(J+)=% .3e  Var(logJ+)=% .3e  step=%.3e\n', ...
            k, Nsteps, Ehist(k), minJhist(k), negHist(k), cvJhist(k), vlogJhist(k), stepHist(k));
    end

    E_opt = Ehist(Nsteps);
end

pUopt = pU;
pVopt = pV;
pBopt = pB;
pTopt = pT;
pLopt = pL;
pRopt = pR;

[Uopt, Vopt] = buildReferenceGridSoftmax(refHeight, pUopt, pVopt, I, J);
BNDopt = buildPhysicalBoundaryFromSoftmax(CURVES, pBopt, pTopt, pLopt, pRopt, I, J);
[Xopt, Yopt, Jc, ~] = solveAnisoLaplaceACFromUV(Uopt, Vopt, BNDopt, A_opt, C_opt, epsDen, hs, ht, derivOrder);

[minq, ESJ] = scaledJacobianStatsHO(Xopt, Yopt, hs, ht, derivOrder, qminSJ);
if strcmpi(filterMode,'minq')
    g_opt = minq;
else
    g_opt = ESJ;
end

info.A_opt = A_opt;
info.C_opt = C_opt;
info.refHeight = refHeight;
info.Ehist = Ehist(1:max(Nsteps,0));
info.minJhist = minJhist(1:max(Nsteps,0));
info.negHist = negHist(1:max(Nsteps,0));
info.cvJhist = cvJhist(1:max(Nsteps,0));
info.vlogJhist = vlogJhist(1:max(Nsteps,0));
info.stepHist = stepHist(1:max(Nsteps,0));
info.minq = minq;
info.ESJ = ESJ;
info.g_opt = g_opt;
info.BNDopt = BNDopt;

fprintf('Final: E=%.6e, minq=%.6e, ESJ=%.6e\n', E_opt, minq, ESJ);

%% -----------------------------------------------------------------------
%% plots
figure; hold on; axis equal; box on;
for jj = 1:J, plot(Xopt(jj,:), Yopt(jj,:), '-k'); end
for ii = 1:I, plot(Xopt(:,ii), Yopt(:,ii), '-k'); end
title(sprintf('Physical grid, simultaneous equalization, A=%.4f, C=%.4f', A_opt, C_opt));
xlabel('X'); ylabel('Y');

figure; hold on; axis equal; box on;
for jj = 1:J, plot(Uopt(jj,:), Vopt(jj,:), '-b'); end
for ii = 1:I, plot(Uopt(:,ii), Vopt(:,ii), '-b'); end
title('Reference grid after simultaneous equalization');
xlabel('u'); ylabel('v');

figure; imagesc(Jc); axis image; colorbar;
title('Physical Jacobian det(dX/dsdt)'); xlabel('i'); ylabel('j');

if Nsteps > 0
    figure; plot(1:Nsteps, info.Ehist, '-o'); grid on;
    xlabel('iter'); ylabel('energy');

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
%% nested objective helpers
    function E = objectiveLogAC(z, pUloc, pVloc, pBloc, pTloc, pLloc, pRloc)
        A = exp(z(1));
        C = exp(z(2));

        if A < A_bracket(1) || A > A_bracket(2) || C < C_bracket(1) || C > C_bracket(2)
            E = 1e50 ...
                + 1e6*(max(0, A_bracket(1)-A)^2 + max(0, A-A_bracket(2))^2) ...
                + 1e6*(max(0, C_bracket(1)-C)^2 + max(0, C-C_bracket(2))^2);
            return;
        end

        E = objectiveAC(A, C, pUloc, pVloc, pBloc, pTloc, pLloc, pRloc);
    end

    function E = objectiveAC(A, C, pUloc, pVloc, pBloc, pTloc, pLloc, pRloc)
        [U, V] = buildReferenceGridSoftmax(refHeight, pUloc, pVloc, I, J);
        BND = buildPhysicalBoundaryFromSoftmax(CURVES, pBloc, pTloc, pLloc, pRloc, I, J);
        [X, Y, Jtmp, ~] = solveAnisoLaplaceACFromUV(U, V, BND, A, C, epsDen, hs, ht, derivOrder);

        [E, ~, ~, ~] = evaluateStateFromJ(X, Y, Jtmp, pUloc, pVloc, pBloc, pTloc, pLloc, pRloc, hs, ht, derivOrder, ...
            thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaRef, betaPhys);
    end
end

%% =======================================================================
function BND = buildPhysicalBoundaryFromSoftmax(CURVES, pB, pT, pL, pR, I, J)
% Reparameterize the four fixed physical boundary curves.

sB = cumulativeSoftmaxSpacing(pB);
sT = cumulativeSoftmaxSpacing(pT);
sL = cumulativeSoftmaxSpacing(pL);
sR = cumulativeSoftmaxSpacing(pR);

[xB, yB] = interpPolylineArc(CURVES.bottomX, CURVES.bottomY, sB);
[xT, yT] = interpPolylineArc(CURVES.topX,    CURVES.topY,    sT);
[xL, yL] = interpPolylineArc(CURVES.leftX,   CURVES.leftY,   sL);
[xR, yR] = interpPolylineArc(CURVES.rightX,  CURVES.rightY,  sR);

BND.X = nan(J,I);
BND.Y = nan(J,I);

BND.X(1,:) = xB(:).';
BND.Y(1,:) = yB(:).';
BND.X(J,:) = xT(:).';
BND.Y(J,:) = yT(:).';
BND.X(:,1) = xL(:);
BND.Y(:,1) = yL(:);
BND.X(:,I) = xR(:);
BND.Y(:,I) = yR(:);

% Enforce exact shared corners.
BND.X(1,1) = CURVES.cBL(1);    BND.Y(1,1) = CURVES.cBL(2);
BND.X(1,I) = CURVES.cBR(1);    BND.Y(1,I) = CURVES.cBR(2);
BND.X(J,1) = CURVES.cTL(1);    BND.Y(J,1) = CURVES.cTL(2);
BND.X(J,I) = CURVES.cTR(1);    BND.Y(J,I) = CURVES.cTR(2);
end

%% =======================================================================
function s = cumulativeSoftmaxSpacing(p)
p = p(:);
d = softmax_safe(p);
d = d / sum(d);
s = [0; cumsum(d)];
s(end) = 1;
end

%% =======================================================================
function q = normalizeLogSpacing(q)
q = q(:);
q = q - mean(q);
end

%% =======================================================================
function [xq, yq] = interpPolylineArc(x, y, tq)
x = x(:);
y = y(:);
tq = tq(:);

seg = sqrt(diff(x).^2 + diff(y).^2);
arc = [0; cumsum(seg)];

if arc(end) <= 0
    xq = repmat(x(1), size(tq));
    yq = repmat(y(1), size(tq));
    return;
end

arc = arc / arc(end);
[arcU, ia] = unique(arc, 'stable');
xU = x(ia);
yU = y(ia);

xq = interp1(arcU, xU, tq, 'linear');
yq = interp1(arcU, yU, tq, 'linear');

xq(1) = x(1);      yq(1) = y(1);
xq(end) = x(end);  yq(end) = y(end);
end

%% =======================================================================
function [scoreUdual, scoreVdual] = logJacobianVarianceScores(Jc, tauJ, eqW, badW)
Jsafe = max(Jc, tauJ);

pos = Jc > 0;
if any(pos(:))
    targetLog = mean(log(max(Jc(pos), tauJ)));
else
    targetLog = log(tauJ);
end

Req = targetLog - log(Jsafe);
Rbad = max(0, tauJ - Jc);

scoreVdual = eqW*mean(Req, 2) + badW*normalizeVec(sum(Rbad, 2));
scoreUdual = eqW*mean(Req, 1).' + badW*normalizeVec(sum(Rbad, 1).');
end

%% =======================================================================
function [E, minJ, negCount, cvPos] = evaluateStateFromJ(X, Y, Jc, pU, pV, pB, pT, pL, pR, hs, ht, derivOrder, ...
    thetaW, mu, lambda, epsJ, kappa, epsJw, gammaOrth, gammaAR, betaRef, betaPhys)

ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder);
EW  = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw);
Ebase = (1-thetaW)*ENH + thetaW*EW;

Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder);
Ear   = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder);

EregRef = betaRef*(sum(pU(:).^2) + sum(pV(:).^2));
EregPhys = betaPhys*(sum(pB(:).^2) + sum(pT(:).^2) + sum(pL(:).^2) + sum(pR(:).^2));

E = Ebase + gammaOrth*Eorth + gammaAR*Ear + EregRef + EregPhys;

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
function tf = isBetterState(newNeg, newMin, newVlog, newE, oldNeg, oldMin, oldVlog, oldE)
tolMin = 1e-14;
tolE = 1e-12;
tf = false;

if newNeg < oldNeg
    tf = true;
elseif newNeg == oldNeg && newMin > oldMin + tolMin
    tf = true;
elseif newNeg == 0 && oldNeg == 0 && newVlog < oldVlog
    tf = true;
elseif newNeg == oldNeg && abs(newMin - oldMin) <= tolMin && newE < oldE - tolE
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

%% =======================================================================
function v = normalizeVec(v)
v = v(:);
nrm = max(abs(v));
if nrm > 0
    v = v / nrm;
end
end

%% =======================================================================
function y = smooth1D(y)
y = y(:);
if numel(y) < 3
    return;
end
ys = y;
ys(2:end-1) = 0.25*y(1:end-2) + 0.5*y(2:end-1) + 0.25*y(3:end);
y = ys;
end

%% =======================================================================
function [Ugrid, Vgrid] = buildReferenceGridSoftmax(refHeight, pU, pV, I, J)
pU = pU(:);
pV = pV(:);

dU = softmax_safe(pU);
dV = softmax_safe(pV);

dU = dU / sum(dU);
dV = dV / sum(dV) * refHeight;

U = [0; cumsum(dU)];
V = [0; cumsum(dV)];

Ugrid = repmat(U(:).', J, 1);
Vgrid = repmat(V(:), 1, I);
end

%% =======================================================================
function y = softmax_safe(x)
x = x(:);
x = x - max(x);
ex = exp(x);
y = ex / sum(ex);
end

%% =======================================================================
function [minq, ESJ] = scaledJacobianStatsHO(X, Y, hs, ht, derivOrder, qmin)
if nargin < 6
    qmin = 0.25;
end

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
function [X, Y, Jc, vlogJ] = solveAnisoLaplaceACFromUV(Ugrid, Vgrid, BND, Acoef, Ccoef, epsDen, hs, ht, derivOrder)
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
Bt = beta_n  / (ht^2);

X = solveOneComponent(BND.X, 1/(Acoef^2));
Y = solveOneComponent(BND.Y, 1/(Ccoef^2));

[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
Jc = xs.*yt - xt.*ys;

pos = Jc > 0;
if any(pos(:))
    vlogJ = var(log(max(Jc(pos), 1e-14)));
else
    vlogJ = inf;
end

    function Z = solveOneComponent(BZ, tScale)
        idx = @(ii,jj) (jj-1)*nI + ii;

        ii = zeros(5*N,1);
        jj = zeros(5*N,1);
        ss = zeros(5*N,1);
        b  = zeros(N,1);
        ptr = 0;

        for jj_int = 1:nJ
            jg = jj_int + 1;
            for ii_int = 1:nI
                ig  = ii_int + 1;
                row = idx(ii_int, jj_int);

                Ww = As(jg, ig-1);
                We = As(jg, ig);
                Ws = tScale * Bt(jg-1, ig);
                Wn = tScale * Bt(jg, ig);

                diagv = Ww + We + Ws + Wn;

                ptr = ptr + 1;
                ii(ptr) = row;
                jj(ptr) = row;
                ss(ptr) = diagv;

                if ii_int > 1
                    ptr = ptr + 1;
                    ii(ptr) = row;
                    jj(ptr) = idx(ii_int-1, jj_int);
                    ss(ptr) = -Ww;
                else
                    b(row) = b(row) + Ww * BZ(jg,1);
                end

                if ii_int < nI
                    ptr = ptr + 1;
                    ii(ptr) = row;
                    jj(ptr) = idx(ii_int+1, jj_int);
                    ss(ptr) = -We;
                else
                    b(row) = b(row) + We * BZ(jg,I);
                end

                if jj_int > 1
                    ptr = ptr + 1;
                    ii(ptr) = row;
                    jj(ptr) = idx(ii_int, jj_int-1);
                    ss(ptr) = -Ws;
                else
                    b(row) = b(row) + Ws * BZ(1,ig);
                end

                if jj_int < nJ
                    ptr = ptr + 1;
                    ii(ptr) = row;
                    jj(ptr) = idx(ii_int, jj_int+1);
                    ss(ptr) = -Wn;
                else
                    b(row) = b(row) + Wn * BZ(J,ig);
                end
            end
        end

        Aop = sparse(ii(1:ptr), jj(1:ptr), ss(1:ptr), N, N);
        Zvec = Aop \ b;

        Z = BZ;
        Z(2:end-1,2:end-1) = reshape(Zvec, [nI, nJ]).';
    end
end

%% =======================================================================
function ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder)
cellArea = hs*ht;
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
Jc = xs.*yt - xt.*ys;
Jsafe = max(Jc, epsJ);
logJ = log(Jsafe);
d = 2;
trC = xs.^2 + ys.^2 + xt.^2 + yt.^2;
W = 0.5*mu*(trC - d) - mu*logJ + 0.5*lambda*(logJ.^2);
r = max(0, epsJ - Jc);
B = kappa*(r.^2);
ENH = cellArea * sum(W(:) + B(:));
end

%% =======================================================================
function EW = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw)
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
a2 = xs.^2 + ys.^2;
b2 = xt.^2 + yt.^2;
Jc = xs.*yt - ys.*xt;
Jsafe = max(Jc, epsJw);
tmp = (a2 + b2) ./ Jsafe;
EW = (hs*ht) * sum(tmp(:));
end

%% =======================================================================
function Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder)
[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
orth = xs.*xt + ys.*yt;
Eorth = (hs*ht) * sum(orth(:).^2);
end

%% =======================================================================
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

%% =======================================================================
function [xs, ys, xt, yt] = metricsHighOrder_nodes(X, Y, hs, ht, derivOrder)
xs = diff1_highorder_dim2(X, hs, derivOrder);
ys = diff1_highorder_dim2(Y, hs, derivOrder);
xt = diff1_highorder_dim1(X, ht, derivOrder);
yt = diff1_highorder_dim1(Y, ht, derivOrder);
end

%% =======================================================================
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

%% =======================================================================
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

%% =======================================================================
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

