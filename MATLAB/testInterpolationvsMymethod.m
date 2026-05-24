%% test_softC_vs_interpolation
% Compare:
%   1. Regular interpolation grid from the same boundary
%   2. Grid Generator output Xopt,Yopt
%
% Reports:
%   E
%   minq
%   ESJ
%   minJ
%   folded cells
%
% IMPORTANT:
%   The Grid Generator function must be saved as:
%
%       gridGenerator.m
%
%   and must be either in this same folder or somewhere on the MATLAB path.

clc; clear; close all;

%% ------------------------------------------------------------------------
% Path setup
thisDir = fileparts(mfilename('fullpath'));

if isempty(thisDir)
    thisDir = pwd;
end

addpath(thisDir);
addpath(fullfile(thisDir,'Grids'));
rehash toolboxcache;

fprintf('Current script folder:\n  %s\n\n', thisDir);

%% ------------------------------------------------------------------------
% Check that MATLAB can find required files/functions

fprintf('Checking required functions...\n');

gridGeneratorPath = which('gridGenerator');
if isempty(gridGeneratorPath)
    error(['MATLAB cannot find gridGenerator.\n\n', ...
           'Fix:\n', ...
           '  1. Save the Grid Generator function in a file named exactly:\n', ...
           '       gridGenerator.m\n', ...
           '  2. Put gridGenerator.m in this folder:\n', ...
           '       %s\n', ...
           '     or add its folder with addpath(...).\n\n', ...
           'Diagnostic command:\n', ...
           '  which gridGenerator -all\n'], thisDir);
else
    fprintf('  Found gridGenerator:\n    %s\n', gridGeneratorPath);
end

makeGridChevronPath = which('makeGridChevron');
if isempty(makeGridChevronPath)
    warning(['MATLAB cannot find makeGridChevron.\n', ...
             'If you use grid_fun = @makeGridChevron, make sure makeGridChevron.m is in the Grids folder or on the path.']);
else
    fprintf('  Found makeGridChevron:\n    %s\n', makeGridChevronPath);
end

makeGridCPath = which('makeGridC');
if ~isempty(makeGridCPath)
    fprintf('  Found makeGridC:\n    %s\n', makeGridCPath);
end

fprintf('Function check complete.\n\n');

%% ------------------------------------------------------------------------
% Parameters

Ns = 41;
Nt = 41;
N_outer = 0010;   % number of ACCEPTED SGD steps with filter

%% ------------------------------------------------------------------------
% Grid function handle

grid_fun = @makeGridC;   % change to @makeGridC if desired

%% ------------------------------------------------------------------------
% Sanity check on grid_fun output sizes

[Xc, Yc] = grid_fun(Ns-1, Nt-1);

fprintf('grid_fun returned size(Xc)=%s, size(Yc)=%s\n', ...
    mat2str(size(Xc)), mat2str(size(Yc)));

if ~isequal(size(Xc), [Nt, Ns]) && isequal(size(Xc), [Ns, Nt])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.';
    Yc = Yc.';
end

if ~isequal(size(Xc), [Nt, Ns])
    error('Expected (Nt x Ns)=(%d x %d), got %s.', ...
        Nt, Ns, mat2str(size(Xc)));
end

%% ------------------------------------------------------------------------
% Run Grid Generator

[C_opt, E_opt, Xopt, Yopt, JcOpt, pUopt, pVopt, Uopt, Vopt, info] = ...
    gridGenerator(grid_fun, Ns, Nt, N_outer);

%% ------------------------------------------------------------------------
% Make sure Grid Generator output has correct orientation

if ~isequal(size(Xopt), [Nt, Ns]) && isequal(size(Xopt), [Ns, Nt])
    warning('Xopt,Yopt returned as (Ns x Nt). Transposing to (Nt x Ns).');
    Xopt = Xopt.';
    Yopt = Yopt.';
end

if ~isequal(size(Xopt), [Nt, Ns])
    error('Expected Xopt,Yopt size (Nt x Ns)=(%d x %d), got %s.', ...
        Nt, Ns, mat2str(size(Xopt)));
end

%% ------------------------------------------------------------------------
% Build regular interpolation grid from the same physical boundary
%
% This uses the boundary of the original grid_fun output.
% The interpolation grid is the ordinary Coons/transfinite interpolation
% baseline, with no optimization.

BND = buildBoundaryFromGrid(Xc, Yc);
[Xinterp, Yinterp] = transfiniteInterpolationGrid(BND);

%% ------------------------------------------------------------------------
% Evaluate both grids

interpStats = evaluateGridStats(Xinterp, Yinterp, Ns, Nt);
genStats    = evaluateGridStats(Xopt,    Yopt,    Ns, Nt);

%% ------------------------------------------------------------------------
% Folded-cell counts

interpFolded = foldedCellMask(interpStats.X, interpStats.Y);
genFolded    = foldedCellMask(genStats.X,    genStats.Y);

nFoldInterp = nnz(interpFolded);
nFoldGen    = nnz(genFolded);

%% ------------------------------------------------------------------------
% Print Grid Generator results

fprintf('\nGRID GENERATOR RESULTS:\n');
fprintf('  C_opt   = %.6f\n', C_opt);
fprintf('  E_opt   = %.6e\n', E_opt);
fprintf('  min(Jc) = %.3e\n', min(JcOpt(:)));

if isfield(info,'accepted')
    fprintf('  accepted steps = %d\n', info.accepted);
end

if isfield(info,'minq')
    fprintf('  minq    = %.6e\n', info.minq);
end

if isfield(info,'ESJ')
    fprintf('  ESJ     = %.6e\n', info.ESJ);
end

%% ------------------------------------------------------------------------
% Print comparison

fprintf('\n============================================================\n');
fprintf('REGULAR INTERPOLATION GRID VS GRID GENERATOR\n');
fprintf('============================================================\n');

fprintf('\nREGULAR INTERPOLATION GRID:\n');
fprintf('  E       = %.6e\n', interpStats.E);
fprintf('  minq    = %.6e\n', interpStats.minq);
fprintf('  ESJ     = %.6e\n', interpStats.ESJ);
fprintf('  minJ    = %.6e\n', min(interpStats.Jc(:)));
fprintf('  folded  = %d\n', nFoldInterp);

fprintf('\nGRID GENERATOR:\n');
fprintf('  E       = %.6e\n', genStats.E);
fprintf('  minq    = %.6e\n', genStats.minq);
fprintf('  ESJ     = %.6e\n', genStats.ESJ);
fprintf('  minJ    = %.6e\n', min(genStats.Jc(:)));
fprintf('  folded  = %d\n', nFoldGen);

fprintf('\nGRID GENERATOR - REGULAR INTERPOLATION:\n');
fprintf('  E change      = %.6e\n', genStats.E - interpStats.E);
fprintf('  minq change   = %.6e\n', genStats.minq - interpStats.minq);
fprintf('  ESJ change    = %.6e\n', genStats.ESJ - interpStats.ESJ);
fprintf('  minJ change   = %.6e\n', min(genStats.Jc(:)) - min(interpStats.Jc(:)));
fprintf('  folded change = %d\n', nFoldGen - nFoldInterp);

fprintf('\nRATIOS, GRID GENERATOR / REGULAR INTERPOLATION:\n');
fprintf('  E ratio       = %.6e\n', genStats.E / max(interpStats.E, eps));
fprintf('  ESJ ratio     = %.6e\n', genStats.ESJ / max(interpStats.ESJ, eps));

%% ------------------------------------------------------------------------
% Plot reference spacings from pUopt,pVopt

dU = exp(pUopt(:));
dU = dU / sum(dU);

dV = exp(pVopt(:));
dV = dV / sum(dV) * C_opt;

figure('Color','w');
subplot(2,1,1);
plot(dU,'-o'); grid on;
title('Reference spacings dU from Grid Generator');
xlabel('i'); ylabel('\Delta U_i');

subplot(2,1,2);
plot(dV,'-o'); grid on;
title(sprintf('Reference spacings dV from Grid Generator, C=%.4f', C_opt));
xlabel('j'); ylabel('\Delta V_j');

%% ------------------------------------------------------------------------
% History plots

if isfield(info,'Ehist') && ~isempty(info.Ehist)
    figure('Color','w');
    plot(info.Ehist,'-o'); grid on;
    title('Grid Generator energy history');
    xlabel('accepted step'); ylabel('E');
end

if isfield(info,'ghist') && ~isempty(info.ghist)
    figure('Color','w');
    plot(info.ghist,'-o'); grid on;
    title('Grid Generator sliver metric history');
    xlabel('accepted step'); ylabel('g');
end

if isfield(info,'minJhist') && ~isempty(info.minJhist)
    figure('Color','w');
    plot(info.minJhist,'-o'); grid on;
    title('Grid Generator min(J) history');
    xlabel('iteration'); ylabel('min(J)');
end

if isfield(info,'vlogJhist') && ~isempty(info.vlogJhist)
    figure('Color','w');
    plot(info.vlogJhist,'-o'); grid on;
    title('Grid Generator Var(log J_+) history');
    xlabel('iteration'); ylabel('Var(log J_+)');
end

%% ------------------------------------------------------------------------
% Plot regular interpolation grid

plotGrid(interpStats.X, interpStats.Y, 'Regular Interpolation Grid');

%% ------------------------------------------------------------------------
% Plot Grid Generator result

plotGrid(genStats.X, genStats.Y, 'Grid Generator');

%% ------------------------------------------------------------------------
% Plot Jacobians side by side

figure('Color','w');

subplot(1,2,1);
imagesc(interpStats.Jc);
axis image; colorbar;
title('Regular Interpolation Grid: Jacobian J');
xlabel('i'); ylabel('j');

subplot(1,2,2);
imagesc(genStats.Jc);
axis image; colorbar;
title('Grid Generator: Jacobian J');
xlabel('i'); ylabel('j');

%% ------------------------------------------------------------------------
% Plot folded cells

plotFoldedCellsComputational(interpStats.X, interpStats.Y, ...
    'Regular Interpolation Grid: folded physical cells');

plotFoldedCellsComputational(genStats.X, genStats.Y, ...
    'Grid Generator: folded physical cells');

%% =======================================================================
function BND = buildBoundaryFromGrid(Xc, Yc)

[J,I] = size(Xc);

BND.X = nan(J,I);
BND.Y = nan(J,I);

BND.X(1,:)   = Xc(1,:);
BND.Y(1,:)   = Yc(1,:);

BND.X(J,:)   = Xc(end,:);
BND.Y(J,:)   = Yc(end,:);

BND.X(:,1)   = Xc(:,1);
BND.Y(:,1)   = Yc(:,1);

BND.X(:,I)   = Xc(:,end);
BND.Y(:,I)   = Yc(:,end);

BND.X(1,1) = Xc(1,1);
BND.Y(1,1) = Yc(1,1);

BND.X(1,I) = Xc(1,end);
BND.Y(1,I) = Yc(1,end);

BND.X(J,1) = Xc(end,1);
BND.Y(J,1) = Yc(end,1);

BND.X(J,I) = Xc(end,end);
BND.Y(J,I) = Yc(end,end);

end

%% =======================================================================
function [X, Y] = transfiniteInterpolationGrid(BND)
% Regular Coons/transfinite interpolation from the four boundary curves.

[J,I] = size(BND.X);

s = linspace(0,1,I);
t = linspace(0,1,J);

[S,T] = meshgrid(s,t);

Xbottom = BND.X(1,:);
Ybottom = BND.Y(1,:);

Xtop = BND.X(end,:);
Ytop = BND.Y(end,:);

Xleft = BND.X(:,1);
Yleft = BND.Y(:,1);

Xright = BND.X(:,end);
Yright = BND.Y(:,end);

X00 = BND.X(1,1);       Y00 = BND.Y(1,1);
X10 = BND.X(1,end);     Y10 = BND.Y(1,end);
X01 = BND.X(end,1);     Y01 = BND.Y(end,1);
X11 = BND.X(end,end);   Y11 = BND.Y(end,end);

Xbottom2 = repmat(Xbottom, J, 1);
Ybottom2 = repmat(Ybottom, J, 1);

Xtop2 = repmat(Xtop, J, 1);
Ytop2 = repmat(Ytop, J, 1);

Xleft2 = repmat(Xleft, 1, I);
Yleft2 = repmat(Yleft, 1, I);

Xright2 = repmat(Xright, 1, I);
Yright2 = repmat(Yright, 1, I);

Xcorner = ...
    (1-S).*(1-T)*X00 + ...
       S .*(1-T)*X10 + ...
    (1-S).*   T *X01 + ...
       S .*   T *X11;

Ycorner = ...
    (1-S).*(1-T)*Y00 + ...
       S .*(1-T)*Y10 + ...
    (1-S).*   T *Y01 + ...
       S .*   T *Y11;

X = ...
    (1-S).*Xleft2 + S.*Xright2 + ...
    (1-T).*Xbottom2 + T.*Xtop2 - ...
    Xcorner;

Y = ...
    (1-S).*Yleft2 + S.*Yright2 + ...
    (1-T).*Ybottom2 + T.*Ytop2 - ...
    Ycorner;

end

%% =======================================================================
function stats = evaluateGridStats(X, Y, Ns, Nt)

I = Ns;
J = Nt;

hs = 1/(I-1);
ht = 1/(J-1);

derivOrder = 2;

thetaW = 0;
mu = 1.0;
lambda = 1.0;
epsJ = 1e-6;
kappa = 10.0;
epsJw = 1e-10;

gammaOrth = 4e-2;
gammaAR   = 4e-2;
betaAB    = 0.0;

qminSJ = 0.25;

u = linspace(0,1,I);
v = linspace(0,1,J);

U = repmat(u, J, 1);
V = repmat(v(:), 1, I);

[xs, ys, xt, yt] = metricsHighOrder(X, Y, hs, ht, derivOrder);
Jc = xs.*yt - xt.*ys;

ENH = neoHookeanEnergyHO(X, Y, hs, ht, mu, lambda, epsJ, kappa, derivOrder);
EW  = winslowEnergyHO(X, Y, hs, ht, derivOrder, epsJw);

Ebase = (1-thetaW)*ENH + thetaW*EW;

Eorth = orthogonalityPenaltyHO(X, Y, hs, ht, derivOrder);
Ear   = aspectRatioPenaltyHO(X, Y, hs, ht, derivOrder);
Ereg  = betaAB * 0.0;

E = Ebase + gammaOrth*Eorth + gammaAR*Ear + Ereg;

[minq, ESJ] = scaledJacobianStatsHO(X, Y, hs, ht, derivOrder, qminSJ);

stats = struct();
stats.E    = E;
stats.minq = minq;
stats.ESJ  = ESJ;

stats.X = X;
stats.Y = Y;
stats.Jc = Jc;
stats.U = U;
stats.V = V;

stats.ENH   = ENH;
stats.EW    = EW;
stats.Eorth = Eorth;
stats.Ear   = Ear;

end

%% =======================================================================
function plotGrid(X, Y, plotTitle)

[J,I] = size(X);

figure('Color','w');
hold on; axis equal; box on;

for j = 1:J
    plot(X(j,:), Y(j,:), '-r');
end

for i = 1:I
    plot(X(:,i), Y(:,i), '-b');
end

title(plotTitle);
xlabel('X'); ylabel('Y');

end

%% =======================================================================
function plotFoldedCellsComputational(X, Y, plotTitle)

[J,I] = size(X);

folded = foldedCellMask(X, Y);

u = linspace(0,1,I);
v = linspace(0,1,J);

figure('Color','w');
hold on; axis equal; box on;
xlim([0 1]); ylim([0 1]);

for j = 1:J-1
    for i = 1:I-1

        px = [u(i), u(i+1), u(i+1), u(i)];
        py = [v(j), v(j), v(j+1), v(j+1)];

        if folded(j,i)
            patch(px, py, 'r', ...
                'FaceAlpha', 0.65, ...
                'EdgeColor', 'k');
        else
            patch(px, py, 'w', ...
                'FaceAlpha', 0.15, ...
                'EdgeColor', [0.75 0.75 0.75]);
        end

    end
end

title(sprintf('%s, folded cells = %d', plotTitle, nnz(folded)));
xlabel('s'); ylabel('t');

end

%% =======================================================================
function folded = foldedCellMask(X, Y)

[J,I] = size(X);
folded = false(J-1,I-1);

tol = 1e-14;

for j = 1:J-1
    for i = 1:I-1

        px = [X(j,i), X(j,i+1), X(j+1,i+1), X(j+1,i)];
        py = [Y(j,i), Y(j,i+1), Y(j+1,i+1), Y(j+1,i)];

        A = 0.5 * sum(px .* py([2 3 4 1]) - py .* px([2 3 4 1]));

        e1 = [px(2)-px(1), py(2)-py(1)];
        e2 = [px(3)-px(2), py(3)-py(2)];
        e3 = [px(4)-px(3), py(4)-py(3)];
        e4 = [px(1)-px(4), py(1)-py(4)];

        c = zeros(4,1);
        c(1) = cross2(e1,e2);
        c(2) = cross2(e2,e3);
        c(3) = cross2(e3,e4);
        c(4) = cross2(e4,e1);

        hasMixedOrientation = any(c > tol) && any(c < -tol);
        hasNonPositiveArea  = A <= tol;

        folded(j,i) = hasMixedOrientation || hasNonPositiveArea;

    end
end

end

%% =======================================================================
function z = cross2(a,b)

z = a(1)*b(2) - a(2)*b(1);

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