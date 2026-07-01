%% test_multigrid_vs_interpolation
% Compare:
%   1. Baseline/interpolation grid returned by grid_fun
%   2. Multigrid Grid Generator output Xopt,Yopt
%
% Reports:
%   E
%   minq
%   ESJ
%   minJ
%   folded cells
%
% IMPORTANT:
%   The multigrid function must be saved as:
%
%       gridGeneratorMultigrid.m
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

gridGeneratorPath = which('gridGeneratorMultigrid');
if isempty(gridGeneratorPath)
    error(['MATLAB cannot find gridGeneratorMultigrid.\n\n', ...
           'Fix:\n', ...
           '  1. Save the Multigrid Grid Generator function in a file named exactly:\n', ...
           '       gridGeneratorMultigrid.m\n', ...
           '  2. Put gridGeneratorMultigrid.m in this folder:\n', ...
           '       %s\n', ...
           '     or add its folder with addpath(...).\n\n', ...
           'Diagnostic command:\n', ...
           '  which gridGeneratorMultigrid -all\n'], thisDir);
else
    fprintf('  Found gridGeneratorMultigrid:\n    %s\n', gridGeneratorPath);
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
N_outer = 10;   % number of multigrid equalization iterations

%% ------------------------------------------------------------------------
% Grid function handle

grid_fun = @makeGridA;   % change to @makeGridC, @makeGridSwan, etc. if desired

%% ------------------------------------------------------------------------
% Get baseline/interpolation grid from grid_fun

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

% Since grid_fun already creates the interpolation/baseline grid,
% use it directly as the comparison grid.
Xinterp = Xc;
Yinterp = Yc;

%% ------------------------------------------------------------------------
% Run Multigrid Grid Generator

[A_opt, C_opt, E_opt, Xopt, Yopt, JcOpt, ...
    pUopt, pVopt, pBopt, pTopt, pLopt, pRopt, ...
    Uopt, Vopt, info] = gridGeneratorMultigrid(grid_fun, Ns, Nt, N_outer);

%% ------------------------------------------------------------------------
% Make sure Multigrid Grid Generator output has correct orientation

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
% Evaluate both grids

interpStats = evaluateGridStats(Xinterp, Yinterp, Ns, Nt);
multiStats  = evaluateGridStats(Xopt,    Yopt,    Ns, Nt);

%% ------------------------------------------------------------------------
% Folded-cell counts

interpFolded = foldedCellMask(interpStats.X, interpStats.Y);
multiFolded  = foldedCellMask(multiStats.X,  multiStats.Y);

nFoldInterp = nnz(interpFolded);
nFoldMulti  = nnz(multiFolded);

%% ------------------------------------------------------------------------
% Print Multigrid Grid Generator results

fprintf('\nMULTIGRID GRID GENERATOR RESULTS:\n');
fprintf('  A_opt   = %.6f\n', A_opt);
fprintf('  C_opt   = %.6f\n', C_opt);
fprintf('  E_opt   = %.6e\n', E_opt);
fprintf('  min(Jc) = %.3e\n', min(JcOpt(:)));

if isfield(info,'minq')
    fprintf('  minq    = %.6e\n', info.minq);
end

if isfield(info,'ESJ')
    fprintf('  ESJ     = %.6e\n', info.ESJ);
end

%% ------------------------------------------------------------------------
% Print comparison

fprintf('\n============================================================\n');
fprintf('BASELINE GRID FROM grid_fun VS MULTIGRID GRID GENERATOR\n');
fprintf('============================================================\n');

fprintf('\nBASELINE GRID FROM grid_fun:\n');
fprintf('  E       = %.6e\n', interpStats.E);
fprintf('  minq    = %.6e\n', interpStats.minq);
fprintf('  ESJ     = %.6e\n', interpStats.ESJ);
fprintf('  minJ    = %.6e\n', min(interpStats.Jc(:)));
fprintf('  folded  = %d\n', nFoldInterp);

fprintf('\nMULTIGRID GRID GENERATOR:\n');
fprintf('  E       = %.6e\n', multiStats.E);
fprintf('  minq    = %.6e\n', multiStats.minq);
fprintf('  ESJ     = %.6e\n', multiStats.ESJ);
fprintf('  minJ    = %.6e\n', min(multiStats.Jc(:)));
fprintf('  folded  = %d\n', nFoldMulti);

fprintf('\nMULTIGRID - BASELINE:\n');
fprintf('  E change      = %.6e\n', multiStats.E - interpStats.E);
fprintf('  minq change   = %.6e\n', multiStats.minq - interpStats.minq);
fprintf('  ESJ change    = %.6e\n', multiStats.ESJ - interpStats.ESJ);
fprintf('  minJ change   = %.6e\n', min(multiStats.Jc(:)) - min(interpStats.Jc(:)));
fprintf('  folded change = %d\n', nFoldMulti - nFoldInterp);

fprintf('\nRATIOS, MULTIGRID / BASELINE:\n');
fprintf('  E ratio       = %.6e\n', multiStats.E / max(interpStats.E, eps));
fprintf('  ESJ ratio     = %.6e\n', multiStats.ESJ / max(interpStats.ESJ, eps));

%% ------------------------------------------------------------------------
% Plot reference spacings from pUopt,pVopt

dU = exp(pUopt(:));
dU = dU / sum(dU);

dV = exp(pVopt(:));
dV = dV / sum(dV);

figure('Color','w');

subplot(2,1,1);
plot(dU,'-o'); grid on;
title('Reference spacing dU from Multigrid Grid Generator');
xlabel('i'); ylabel('\Delta U_i');

subplot(2,1,2);
plot(dV,'-o'); grid on;
title('Reference spacing dV from Multigrid Grid Generator');
xlabel('j'); ylabel('\Delta V_j');

%% ------------------------------------------------------------------------
% Plot physical boundary spacings from pBopt,pTopt,pLopt,pRopt

dB = exp(pBopt(:)); dB = dB / sum(dB);
dT = exp(pTopt(:)); dT = dT / sum(dT);
dL = exp(pLopt(:)); dL = dL / sum(dL);
dR = exp(pRopt(:)); dR = dR / sum(dR);

figure('Color','w');

subplot(2,2,1);
plot(dB,'-o'); grid on;
title('Bottom boundary spacing dB');
xlabel('i'); ylabel('\Delta s^B_i');

subplot(2,2,2);
plot(dT,'-o'); grid on;
title('Top boundary spacing dT');
xlabel('i'); ylabel('\Delta s^T_i');

subplot(2,2,3);
plot(dL,'-o'); grid on;
title('Left boundary spacing dL');
xlabel('j'); ylabel('\Delta s^L_j');

subplot(2,2,4);
plot(dR,'-o'); grid on;
title('Right boundary spacing dR');
xlabel('j'); ylabel('\Delta s^R_j');

%% ------------------------------------------------------------------------
% History plots

if isfield(info,'Ehist') && ~isempty(info.Ehist)
    figure('Color','w');
    plot(info.Ehist,'-o'); grid on;
    title('Multigrid Grid Generator energy history');
    xlabel('iteration'); ylabel('E');
end

if isfield(info,'minJhist') && ~isempty(info.minJhist)
    figure('Color','w');
    plot(info.minJhist,'-o'); grid on;
    title('Multigrid Grid Generator min(J) history');
    xlabel('iteration'); ylabel('min(J)');
end

if isfield(info,'negHist') && ~isempty(info.negHist)
    figure('Color','w');
    plot(info.negHist,'-o'); grid on;
    title('Multigrid Grid Generator invalid-cell history');
    xlabel('iteration'); ylabel('# cells with J <= 0');
end

if isfield(info,'cvJhist') && ~isempty(info.cvJhist)
    figure('Color','w');
    plot(info.cvJhist,'-o'); grid on;
    title('Multigrid Grid Generator CV(J_+) history');
    xlabel('iteration'); ylabel('CV(J_+)');
end

if isfield(info,'vlogJhist') && ~isempty(info.vlogJhist)
    figure('Color','w');
    plot(info.vlogJhist,'-o'); grid on;
    title('Multigrid Grid Generator Var(log J_+) history');
    xlabel('iteration'); ylabel('Var(log J_+)');
end

if isfield(info,'stepHist') && ~isempty(info.stepHist)
    figure('Color','w');
    plot(info.stepHist,'-o'); grid on;
    title('Multigrid Grid Generator accepted step history');
    xlabel('iteration'); ylabel('accepted step size');
end

%% ------------------------------------------------------------------------
% Plot baseline grid

plotGrid(interpStats.X, interpStats.Y, 'Baseline Grid from grid\_fun');

%% ------------------------------------------------------------------------
% Plot Multigrid Grid Generator result

plotGrid(multiStats.X, multiStats.Y, ...
    sprintf('Multigrid Grid Generator, A=%.4f, C=%.4f', A_opt, C_opt));

%% ------------------------------------------------------------------------
% Plot reference grid returned by multigrid method

plotGrid(Uopt, Vopt, 'Final Reference Grid from Multigrid Generator');

%% ------------------------------------------------------------------------
% Plot Jacobians side by side

figure('Color','w');

subplot(1,2,1);
imagesc(interpStats.Jc);
axis image; colorbar;
title('Baseline Grid: Jacobian J');
xlabel('i'); ylabel('j');

subplot(1,2,2);
imagesc(multiStats.Jc);
axis image; colorbar;
title('Multigrid Grid Generator: Jacobian J');
xlabel('i'); ylabel('j');

%% ------------------------------------------------------------------------
% Plot folded cells

plotFoldedCellsComputational(interpStats.X, interpStats.Y, ...
    'Baseline Grid: folded physical cells');

plotFoldedCellsComputational(multiStats.X, multiStats.Y, ...
    'Multigrid Grid Generator: folded physical cells');