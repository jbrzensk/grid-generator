%% mole_mg_fish_tfi_first.m
% Directly runnable script.
%
% Compare MOLE curvilinear heat solves on:
%
%   1. Transfinite interpolation grid built from the Fish-grid boundary
%   2. Multigrid gridGeneratorMultigrid output built from the same Fish grid
%
% This is intentionally a SCRIPT, not a function. Put this file in
% mole3/grid-generator/MATLAB and press Run.
%
% No folded-grid / signed-cell-area checks are performed. The script simply
% tries the transfinite-interpolation solve first, then the multigrid-generator solve.
% If a solve errors, continueOnSolverError controls whether the script records
% the failure in the table or stops immediately.

clc; clear; close all;

%% ------------------------------------------------------------------------
% Path setup

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end

addpath(thisDir);
addpath(fullfile(thisDir, 'Grids'));

% Common MOLE path layouts. Nonexistent folders are harmless on the path.
addpath(fullfile(thisDir, 'src', 'matlab_octave'));
addpath(fullfile(thisDir, '..', 'src', 'matlab_octave'));
addpath(fullfile(thisDir, '..', '..', 'src', 'matlab_octave'));
addpath('../../src/matlab_octave');

rehash toolboxcache;

%% ------------------------------------------------------------------------
% Parameters you are most likely to edit

% Fish input grid. This local no-plot version mirrors Grids/makeGridFish.m,
% but avoids timing/plot clutter from the grid input script.
grid_fun = @makeGridFishNoPlot;
gridName = 'Fish';

% Multigrid generator is usually easiest to view on a square grid.
Ns = 11;
Nt = 11;

% One rectangular run:
gridSizes = [Ns, Nt];

% For several runs, use one [Ns Nt] pair per row, for example:
% gridSizes = [21 21; 41 41; 81 81];
% gridSizes = [41 21; 81 41; 121 61];

plotGridSize = [Ns, Nt];      % [] disables solution/grid plots; [Ns Nt] plots that one run
Nsteps = 50;                  % gridGeneratorMultigrid SGD/equalization steps

%% ------------------------------------------------------------------------
% MOLE / PDE parameters

k = 2;                        % MOLE order of accuracy
alpha = 1.0;                  % diffusion coefficient
method = 'implicit';          % 'implicit' or 'explicit'

Tfinal = 1.0;
CFL = 0.2;
plotEvery = 25;

%% ------------------------------------------------------------------------
% Manufactured solution for arbitrary physical coordinates
%
% u(x,y,t) = exp(-t)*cos(exactA*x)*sin(exactB*y)
% sourceTerm is chosen so this is the exact solution of
% u_t = alpha*Lap(u) + f.

exactA = pi;
exactB = pi;

%% ------------------------------------------------------------------------
% Plot / output options

makePlots = true;
plotGrids = true;
suppressGridFunctionFigures = true;

% This is not a folded-grid check. It just lets the generated-grid result
% and the output table survive if either solve errors in MOLE.
continueOnSolverError = true;

saveResults = true;
outputDir = pwd;
outputPrefix = '';            % leave empty to auto-name from gridName

%% ------------------------------------------------------------------------
% Check required functions

if isempty(which('gridGeneratorMultigrid'))
    error(['Cannot find gridGeneratorMultigrid.m. Put it on the MATLAB path or in this folder.\n', ...
           'Diagnostic: run "which gridGeneratorMultigrid -all".']);
end

if isempty(which('nodal2DCurv'))
    error(['Cannot find nodal2DCurv.m. Make sure MOLE src/matlab_octave is on the path.\n', ...
           'Diagnostic: run "which nodal2DCurv -all".']);
end

if isempty(grid_fun) || ~isa(grid_fun, 'function_handle')
    error('grid_fun must be a function handle, for example grid_fun = @makeGridFishNoPlot;');
end

method = validatestringCompat(method, {'explicit', 'implicit'}, 'method');
gridSizes = normalizeGridSizes(gridSizes);

gridDescriptor = func2str(grid_fun);
if isempty(gridName)
    gridName = gridDescriptor;
else
    gridName = char(gridName);
end

if isempty(outputPrefix)
    outputPrefix = ['mole_' sanitizeFileStem(gridName) '_TFIfirst_multigridGridGenerator_vs_transfinite_results'];
else
    outputPrefix = char(outputPrefix);
end

fprintf('Using input grid handle:\n  %s\n', gridDescriptor);
fprintf('Using gridGeneratorMultigrid:\n  %s\n', which('gridGeneratorMultigrid'));
fprintf('Using nodal2DCurv:\n  %s\n\n', which('nodal2DCurv'));

%% ------------------------------------------------------------------------
% Storage for table/results

rows = {};
summaryRows = {};
RunResults = struct([]);

%% ------------------------------------------------------------------------
% Loop over grid sizes

for rr = 1:size(gridSizes, 1)

    Ns = gridSizes(rr, 1);
    Nt = gridSizes(rr, 2);

    makePlotsThisRun = shouldMakePlots(Ns, Nt, plotGridSize, makePlots);

    fprintf('\n============================================================\n');
    fprintf('RUNNING %s GRID SIZE %d x %d\n', gridName, Ns, Nt);
    fprintf('============================================================\n');

    %% --------------------------------------------------------------------
    % Build initial boundary-fitted grid from the Fish grid function.
    % The grid generator wants a function taking cell counts, so use
    % Ns-1,Nt-1 here.

    [Xc, Yc] = callGridFunction(grid_fun, Ns-1, Nt-1, suppressGridFunctionFigures);

    fprintf('grid_fun returned size(Xc)=%s, size(Yc)=%s\n', ...
        mat2str(size(Xc)), mat2str(size(Yc)));

    [Xc, Yc] = normalizeGridOrientation(Xc, Yc, Nt, Ns, 'input grid function');

    %% --------------------------------------------------------------------
    % Build transfinite interpolation baseline from the input-grid boundary.

    BND = buildBoundaryFromGrid(Xc, Yc);
    [Xtfi, Ytfi] = transfiniteInterpolationGrid(BND);

    if makePlotsThisRun && plotGrids
        plotGridLines(Xtfi, Ytfi, sprintf('%s transfinite interpolation grid, %d by %d', gridName, Ns, Nt));
    end

    %% --------------------------------------------------------------------
    % Solve PDE on transfinite interpolation grid first.
    % There are no folded-grid checks before this call.

    resTFI = runMoleHeatSolve( ...
        Xtfi, Ytfi, ...
        'Transfinite interpolation', ...
        k, alpha, method, Tfinal, CFL, plotEvery, ...
        exactA, exactB, makePlotsThisRun, continueOnSolverError);

    %% --------------------------------------------------------------------
    % Generate optimized grid using multigrid gridGeneratorMultigrid.
    % Suppress all figures created internally by gridGeneratorMultigrid.

    oldFigVis = get(0, 'DefaultFigureVisible');
    figsBefore = findall(0, 'Type', 'figure');
    set(0, 'DefaultFigureVisible', 'off');

    try
        [A_opt, C_opt, E_grid, Xgen, Ygen, Jc, ...
            pUopt, pVopt, pBopt, pTopt, pLopt, pRopt, ...
            Uref, Vref, gridInfo] = gridGeneratorMultigrid(grid_fun, Ns, Nt, Nsteps);
    catch ME
        set(0, 'DefaultFigureVisible', oldFigVis);
        deleteNewFigures(figsBefore);
        rethrow(ME);
    end

    set(0, 'DefaultFigureVisible', oldFigVis);
    deleteNewFigures(figsBefore);

    [Xgen, Ygen] = normalizeGridOrientation(Xgen, Ygen, Nt, Ns, 'gridGeneratorMultigrid');

    fprintf('\nGenerated multigrid grid complete for %d x %d:\n', Ns, Nt);
    fprintf('  A_opt       = %.6e\n', A_opt);
    fprintf('  C_opt       = %.6e\n', C_opt);
    fprintf('  Grid energy = %.6e\n', E_grid);
    fprintf('  min(Jc)     = %.6e\n', min(Jc(:)));
    fprintf('  max(Jc)     = %.6e\n', max(Jc(:)));

    if isfield(gridInfo, 'accepted')
        fprintf('  accepted    = %d\n', gridInfo.accepted);
    end
    if isfield(gridInfo, 'thetaW')
        fprintf('  thetaW      = %.6e\n', gridInfo.thetaW);
    end
    if isfield(gridInfo, 'minq')
        fprintf('  minq        = %.6e\n', gridInfo.minq);
    end
    if isfield(gridInfo, 'ESJ')
        fprintf('  ESJ         = %.6e\n', gridInfo.ESJ);
    end

    if makePlotsThisRun && plotGrids
        plotGridLines(Xgen, Ygen, sprintf('%s multigrid gridGenerator grid, %d by %d', gridName, Ns, Nt));
    end

    %% --------------------------------------------------------------------
    % Compare multigrid gridGenerator geometry against the transfinite grid.
    % Note: the multigrid method may reparameterize the physical boundaries,
    % so this is a node-to-node grid difference, not a boundary-curve error.

    dGrid = hypot(Xgen - Xtfi, Ygen - Ytfi);
    rmsDistToTFI = sqrt(mean(dGrid(:).^2));
    maxDistToTFI = max(dGrid(:));

    fprintf('\nGeometry difference: gridGeneratorMultigrid vs transfinite interpolation\n');
    fprintf('  RMS distance = %.6e\n', rmsDistToTFI);
    fprintf('  Max distance = %.6e\n', maxDistToTFI);

    %% --------------------------------------------------------------------
    % Solve PDE on multigrid gridGenerator grid second.

    resGEN = runMoleHeatSolve( ...
        Xgen, Ygen, ...
        'Multigrid gridGenerator', ...
        k, alpha, method, Tfinal, CFL, plotEvery, ...
        exactA, exactB, makePlotsThisRun, continueOnSolverError);

    %% --------------------------------------------------------------------
    % Save rows for table.

    gridSizeText = sprintf('%d x %d', Ns, Nt);

    rows(end+1,:) = makeRow(gridName, gridSizeText, 'Transfinite interpolation', ...
        resTFI, NaN, NaN, NaN, 0, 0); %#ok<SAGROW>

    rows(end+1,:) = makeRow(gridName, gridSizeText, 'Multigrid gridGenerator', ...
        resGEN, A_opt, C_opt, E_grid, rmsDistToTFI, maxDistToTFI); %#ok<SAGROW>

    summaryRows(end+1,:) = makeSummaryRow(gridName, gridSizeText, resGEN, resTFI); %#ok<SAGROW>

    RunResults(rr).gridName = gridName; %#ok<SAGROW>
    RunResults(rr).Ns = Ns;
    RunResults(rr).Nt = Nt;
    RunResults(rr).inputGrid.X = Xc;
    RunResults(rr).inputGrid.Y = Yc;
    RunResults(rr).transfiniteGrid.X = Xtfi;
    RunResults(rr).transfiniteGrid.Y = Ytfi;
    RunResults(rr).generatedGrid.X = Xgen;
    RunResults(rr).generatedGrid.Y = Ygen;
    RunResults(rr).generatedGrid.Aopt = A_opt;
    RunResults(rr).generatedGrid.Copt = C_opt;
    RunResults(rr).generatedGrid.Egrid = E_grid;
    RunResults(rr).generatedGrid.Jc = Jc;
    RunResults(rr).generatedGrid.pUopt = pUopt;
    RunResults(rr).generatedGrid.pVopt = pVopt;
    RunResults(rr).generatedGrid.pBopt = pBopt;
    RunResults(rr).generatedGrid.pTopt = pTopt;
    RunResults(rr).generatedGrid.pLopt = pLopt;
    RunResults(rr).generatedGrid.pRopt = pRopt;
    RunResults(rr).generatedGrid.Uref = Uref;
    RunResults(rr).generatedGrid.Vref = Vref;
    RunResults(rr).generatedGrid.info = gridInfo;
    RunResults(rr).resGEN = resGEN;
    RunResults(rr).resTFI = resTFI;
    RunResults(rr).rmsDistToTFI = rmsDistToTFI;
    RunResults(rr).maxDistToTFI = maxDistToTFI;

end

%% ------------------------------------------------------------------------
% Print full table

fprintf('\n\n============================================================\n');
fprintf('FULL MOLE CURVILINEAR HEAT SOLVE TABLE\n');
fprintf('============================================================\n');

fprintf('\n%-12s %-12s %-28s %-10s %14s %10s %14s %14s %14s\n', ...
    'InputGrid', 'Grid', 'Method', 'Status', 'dt', 'Steps', 'L2', 'Linf', 'diagRate');

fprintf('%s\n', repmat('-', 1, 128));

for r = 1:size(rows,1)
    fprintf('%-12s %-12s %-28s %-10s %14.6e %10d %14.6e %14.6e %14.6e\n', ...
        rows{r,1}, rows{r,2}, rows{r,3}, rows{r,4}, ...
        rows{r,5}, rows{r,6}, rows{r,7}, rows{r,8}, rows{r,9});
end

fprintf('\n\n============================================================\n');
fprintf('IMPROVEMENT SUMMARY\n');
fprintf('============================================================\n');

fprintf('\n%-12s %-12s %14s %14s  %s\n', ...
    'InputGrid', 'Grid', 'L2_TFI/L2_GEN', 'Linf_TFI/Linf_GEN', 'Conclusion');

fprintf('%s\n', repmat('-', 1, 100));

for r = 1:size(summaryRows,1)
    fprintf('%-12s %-12s %14.6e %14.6e  %s\n', ...
        summaryRows{r,1}, summaryRows{r,2}, summaryRows{r,3}, summaryRows{r,4}, summaryRows{r,5});
end

%% ------------------------------------------------------------------------
% MATLAB tables and saved files

ResultsTable = cell2table(rows, ...
    'VariableNames', {'InputGrid','Grid','Method','Status','dt','Steps','L2','Linf', ...
                      'diagRate','Aopt','Copt','GridEnergy','RMSDistToTFI','MaxDistToTFI','Message'});

SummaryTable = cell2table(summaryRows, ...
    'VariableNames', {'InputGrid','Grid','L2_TFI_over_L2_GEN','Linf_TFI_over_Linf_GEN','Conclusion'});

disp(ResultsTable);
disp(SummaryTable);

if saveResults
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    csvFile = fullfile(outputDir, [outputPrefix '.csv']);
    xlsxFile = fullfile(outputDir, [outputPrefix '.xlsx']);
    summaryCsvFile = fullfile(outputDir, [outputPrefix '_summary.csv']);
    summaryXlsxFile = fullfile(outputDir, [outputPrefix '_summary.xlsx']);
    matFile = fullfile(outputDir, [outputPrefix '.mat']);

    writetable(ResultsTable, csvFile);
    writetable(SummaryTable, summaryCsvFile);

    wroteXlsx = false;
    try
        writetable(ResultsTable, xlsxFile);
        writetable(SummaryTable, summaryXlsxFile);
        wroteXlsx = true;
    catch ME
        warning('MOLECompare:WriteXlsxFailed', 'Could not write xlsx files: %s', ME.message);
    end

    save(matFile, 'ResultsTable', 'SummaryTable', 'RunResults');

    fprintf('\nSaved results to:\n');
    fprintf('  %s\n', csvFile);
    if wroteXlsx
        fprintf('  %s\n', xlsxFile);
    end
    fprintf('  %s\n', summaryCsvFile);
    if wroteXlsx
        fprintf('  %s\n', summaryXlsxFile);
    end
    fprintf('  %s\n', matFile);
end

%% ========================================================================
function row = makeRow(inputGridName, gridSizeText, methodName, res, A_opt, C_opt, E_grid, rmsDistToTFI, maxDistToTFI)

row = {char(inputGridName), char(gridSizeText), char(methodName), char(res.Status), ...
       res.dt, res.ntSteps, res.L2err, res.Linferr, res.diagRate, ...
       A_opt, C_opt, E_grid, rmsDistToTFI, maxDistToTFI, char(res.Message)};

end

%% ========================================================================
function row = makeSummaryRow(inputGridName, gridSizeText, resGEN, resTFI)

if strcmpi(resGEN.Status, 'ok') && strcmpi(resTFI.Status, 'ok')
    L2ratio = resTFI.L2err / resGEN.L2err;
    Linfratio = resTFI.Linferr / resGEN.Linferr;

    if L2ratio > 1 && Linfratio > 1
        conclusion = sprintf('multigrid gridGenerator is better: L2 %.3gx lower, Linf %.3gx lower', L2ratio, Linfratio);
    elseif L2ratio > 1
        conclusion = sprintf('multigrid gridGenerator has lower L2 by %.3gx', L2ratio);
    elseif Linfratio > 1
        conclusion = sprintf('multigrid gridGenerator has lower Linf by %.3gx', Linfratio);
    else
        conclusion = 'transfinite interpolation has equal/lower reported error';
    end

elseif strcmpi(resGEN.Status, 'ok') && ~strcmpi(resTFI.Status, 'ok')
    L2ratio = NaN;
    Linfratio = NaN;
    conclusion = ['multigrid gridGenerator solved; transfinite interpolation failed: ' char(resTFI.Message)];

elseif ~strcmpi(resGEN.Status, 'ok') && strcmpi(resTFI.Status, 'ok')
    L2ratio = NaN;
    Linfratio = NaN;
    conclusion = ['transfinite interpolation solved; multigrid gridGenerator failed: ' char(resGEN.Message)];

else
    L2ratio = NaN;
    Linfratio = NaN;
    conclusion = 'both solves failed';
end

row = {char(inputGridName), char(gridSizeText), L2ratio, Linfratio, conclusion};

end

%% ========================================================================
function result = runMoleHeatSolve(X, Y, gridLabel, k, alpha, method, Tfinal, CFL, plotEvery, exactA, exactB, makePlots, continueOnSolverError)

try
    result = solveMoleHeatOnGrid(X, Y, gridLabel, k, alpha, method, Tfinal, CFL, plotEvery, exactA, exactB, makePlots);
    result.Status = 'ok';
    result.Message = '';
catch ME
    if ~continueOnSolverError
        rethrow(ME);
    end

    warning('MOLECompare:SolveFailed', 'MOLE solve failed on %s: %s', gridLabel, ME.message);

    result = failedResult(gridLabel, X, Y, ME);
end

end

%% ========================================================================
function result = failedResult(gridLabel, X, Y, ME)

[Nt, Ns] = size(X);

result = struct();
result.gridName = char(gridLabel);
result.Status = 'failed';
result.Message = ME.message;
result.dt = NaN;
result.ntSteps = NaN;
result.L2err = NaN;
result.Linferr = NaN;
result.X = X;
result.Y = Y;
result.U = NaN(size(X));
result.Uexact = NaN(size(X));
result.Err = NaN(size(X));
result.diagRate = NaN;
result.Ns = Ns;
result.Nt = Nt;

end

%% ========================================================================
function result = solveMoleHeatOnGrid(X, Y, gridLabel, k, alpha, method, Tfinal, CFL, plotEvery, exactA, exactB, makePlots)

[Nt, Ns] = size(X);
[n, m] = size(X);
Nnodes = n*m;

%% Build MOLE curvilinear nodal Laplacian

[Nx, Ny] = nodal2DCurv(k, X, Y);
L = [Nx Ny]*[Nx; Ny];

%% Boundary indices

bdry_idx = boundaryIdx2D(m, n);

%% CFL-style timestep

diagRateVec = full(abs(diag(L)));
diagRate = alpha * max(diagRateVec(:));

if diagRate <= 0 || ~isfinite(diagRate)
    error('Invalid CFL rate for %s. Check the curvilinear operator L.', gridLabel);
end

dt = CFL / diagRate;
ntSteps = ceil(Tfinal / dt);
dt = Tfinal / ntSteps;

fprintf('\nMOLE curvilinear heat solve on %s, %d x %d\n', gridLabel, Ns, Nt);
fprintf('  method     = %s\n', method);
fprintf('  nodes      = %d\n', Nnodes);
fprintf('  grid       = %d x %d\n', Ns, Nt);
fprintf('  max rate   = %.6e\n', diagRate);
fprintf('  CFL        = %.6e\n', CFL);
fprintf('  dt         = %.6e\n', dt);
fprintf('  time steps = %d\n', ntSteps);

%% Initial condition

t = 0.0;

Umat0 = uExact(X, Y, t, exactA, exactB);
Umat  = Umat0;

U = matToVec(Umat);

BCmat = uExact(X, Y, t, exactA, exactB);
BCvec = matToVec(BCmat);
U(bdry_idx) = BCvec(bdry_idx);

%% Time-stepping matrix

Iop = speye(Nnodes);

switch lower(method)
    case 'implicit'
        A = Iop - dt*alpha*L;

        A(bdry_idx,:) = 0;
        for q = 1:numel(bdry_idx)
            A(bdry_idx(q), bdry_idx(q)) = 1;
        end

    case 'explicit'
        A = Iop + dt*alpha*L;

    otherwise
        error('method must be "implicit" or "explicit".');
end

%% Time integration

for step = 1:ntSteps

    switch lower(method)

        case 'implicit'
            tnew = t + dt;

            Fmat = sourceTerm(X, Y, tnew, exactA, exactB, alpha);
            rhs = U + dt*matToVec(Fmat);

            BCmat = uExact(X, Y, tnew, exactA, exactB);
            BCvec = matToVec(BCmat);
            rhs(bdry_idx) = BCvec(bdry_idx);

            U = A\rhs;
            t = tnew;

        case 'explicit'
            Fmat = sourceTerm(X, Y, t, exactA, exactB, alpha);
            rhs = matToVec(Fmat);

            U = A*U + dt*rhs;
            t = t + dt;

            BCmat = uExact(X, Y, t, exactA, exactB);
            BCvec = matToVec(BCmat);
            U(bdry_idx) = BCvec(bdry_idx);
    end

    if mod(step,plotEvery) == 0 || step == ntSteps
        fprintf('  step %5d / %5d, t = %.5f\n', step, ntSteps, t);
    end
end

Umat = vecToMat(U, n, m);

%% Exact solution and error

Uexact = uExact(X, Y, Tfinal, exactA, exactB);
Err = Umat - Uexact;

L2err = sqrt(mean(Err(:).^2));
Linferr = max(abs(Err(:)));

fprintf('\nERROR for %s, %d x %d at T = %.5f:\n', gridLabel, Ns, Nt, Tfinal);
fprintf('  L2 nodal   = %.6e\n', L2err);
fprintf('  Linf nodal = %.6e\n', Linferr);

%% Shared color limits

solutionValues = [Umat0(:); Umat(:); Uexact(:)];
solCLim = [min(solutionValues), max(solutionValues)];

if abs(solCLim(2) - solCLim(1)) < 1e-14
    solCLim = solCLim + [-1, 1]*1e-14;
end

%% Plots only for selected grid size

if makePlots

    safeLabel = char(gridLabel);

    plotPhysicalHeatmap(X, Y, Umat0, ...
        sprintf('%s, %d by %d: initial condition', safeLabel, Ns, Nt), ...
        'temperature', solCLim);

    plotPhysicalHeatmap(X, Y, Umat, ...
        sprintf('Numerical solution on %s at t=%.1f', safeLabel, Tfinal), ...
        'temperature', solCLim);

    plotPhysicalHeatmap(X, Y, Uexact, ...
        sprintf('Exact solution on %s at t=%.1f', safeLabel, Tfinal), ...
        'temperature', solCLim);

    plotPhysicalHeatmap(X, Y, Err, ...
        sprintf('Error on %s at t=%.1f: L_2 = %.3e, L_\infty = %.3e', safeLabel, Tfinal, L2err, Linferr), ...
        'error', []);
end

%% Return result

result = struct();
result.gridName = char(gridLabel);
result.dt = dt;
result.ntSteps = ntSteps;
result.L2err = L2err;
result.Linferr = Linferr;
result.X = X;
result.Y = Y;
result.U = Umat;
result.Uexact = Uexact;
result.Err = Err;
result.diagRate = diagRate;
result.Ns = Ns;
result.Nt = Nt;

end

%% ========================================================================
function BND = buildBoundaryFromGrid(Xc, Yc)

[J,I] = size(Xc);

BND.X = nan(J,I);
BND.Y = nan(J,I);

BND.X(1,:) = Xc(1,:);
BND.Y(1,:) = Yc(1,:);

BND.X(J,:) = Xc(end,:);
BND.Y(J,:) = Yc(end,:);

BND.X(:,1) = Xc(:,1);
BND.Y(:,1) = Yc(:,1);

BND.X(:,I) = Xc(:,end);
BND.Y(:,I) = Yc(:,end);

BND.X(1,1) = Xc(1,1);
BND.Y(1,1) = Yc(1,1);

BND.X(1,I) = Xc(1,end);
BND.Y(1,I) = Yc(1,end);

BND.X(J,1) = Xc(end,1);
BND.Y(J,1) = Yc(end,1);

BND.X(J,I) = Xc(end,end);
BND.Y(J,I) = Yc(end,end);

end

%% ========================================================================
function [X, Y] = transfiniteInterpolationGrid(BND)

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

%% ========================================================================
function [X, Y] = makeGridFishNoPlot(nptsx, nptsy)
% No-plot version of Grids/makeGridFish.m.

Ni = nptsx;
Nj = nptsy;

xi  = linspace(0,1,Ni+1);
eta = linspace(0,1,Nj+1);

[X, Y] = meshgrid(xi, eta);

X(2:end-1, 2:end-1) = 0;
Y(2:end-1, 2:end-1) = 0;

%% Right bulge inward
xpt = X(:,end);
ypt = Y(:,end);

amp = -0.1;
xnew = amp .* cos(2*pi .* ypt) - amp;
xnew = xpt - xnew;
X(:,end) = xnew;

%% Bottom dome
xpt = X(1,:);
amp = -0.2;
ynew = amp .* cos(2*pi .* xpt) - amp;
Y(1,:) = ynew;

%% Top valley
xpt = X(end,:);
amp = 0.2;
ynew = (amp .* cos(2*pi .* xpt)) + (1 - amp);
Y(end,:) = ynew;

%% Quasi-interpolation interior
for j = 1:Nj+1
    t = eta(j);
    for i = 1:Ni+1
        s = xi(i);
        X(j,i) = (1 - s) * X(j,1) + s * X(j,end);
        Y(j,i) = (1 - t) * Y(1,i) + t * Y(end,i);
    end
end

end

%% ========================================================================
function u = uExact(x, y, t, exactA, exactB)

u = exp(-t) .* cos(exactA*x) .* sin(exactB*y);

end

%% ========================================================================
function f = sourceTerm(x, y, t, exactA, exactB, alpha)

u = uExact(x, y, t, exactA, exactB);
f = (alpha*(exactA^2 + exactB^2) - 1) .* u;

end

%% ========================================================================
function v = matToVec(Umat)

v = reshape(Umat.', [], 1);

end

%% ========================================================================
function Umat = vecToMat(v, n, m)

Umat = reshape(v, m, n).';

end

%% ========================================================================
function idx = boundaryIdx2D(m, n)

idx = [];

for i = 1:m
    idx(end+1,1) = sub2vec(1, i, m); %#ok<AGROW>
    idx(end+1,1) = sub2vec(n, i, m); %#ok<AGROW>
end

for j = 2:n-1
    idx(end+1,1) = sub2vec(j, 1, m); %#ok<AGROW>
    idx(end+1,1) = sub2vec(j, m, m); %#ok<AGROW>
end

idx = unique(idx);

end

%% ========================================================================
function k = sub2vec(j, i, m)

k = (j-1)*m + i;

end

%% ========================================================================
function plotPhysicalHeatmap(X, Y, Q, plotTitle, colorbarLabel, climVals)

figure('Color','w');

surf(X, Y, Q, 'EdgeColor', 'none');
view(2);
axis equal tight;
box on;

cb = colorbar;

if nargin >= 6 && ~isempty(climVals)
    applyCLim(climVals);
end

ylabel(cb, colorbarLabel);
title(plotTitle, 'Interpreter', 'none');
xlabel('x');
ylabel('y');

end

%% ========================================================================
function plotGridLines(X, Y, plotTitle)

figure('Color','w');
hold on;

for j = 1:size(X,1)
    plot(X(j,:), Y(j,:), 'k-', 'LineWidth', 0.8);
end

for i = 1:size(X,2)
    plot(X(:,i), Y(:,i), 'k-', 'LineWidth', 0.8);
end

axis equal tight;
box on;
title(plotTitle, 'Interpreter', 'none');
xlabel('x');
ylabel('y');

end

%% ========================================================================
function [X, Y] = callGridFunction(grid_fun, nptsx, nptsy, suppressFigures)

if ~suppressFigures
    [X, Y] = grid_fun(nptsx, nptsy);
    return;
end

oldFigVis = get(0, 'DefaultFigureVisible');
figsBefore = findall(0, 'Type', 'figure');
set(0, 'DefaultFigureVisible', 'off');

try
    [X, Y] = grid_fun(nptsx, nptsy);
catch ME
    set(0, 'DefaultFigureVisible', oldFigVis);
    deleteNewFigures(figsBefore);
    rethrow(ME);
end

set(0, 'DefaultFigureVisible', oldFigVis);
deleteNewFigures(figsBefore);

end

%% ========================================================================
function deleteNewFigures(figsBefore)

figsAfter = findall(0, 'Type', 'figure');
newFigs = setdiff(figsAfter, figsBefore);
if ~isempty(newFigs)
    delete(newFigs);
end

end

%% ========================================================================
function [X, Y] = normalizeGridOrientation(X, Y, Nt, Ns, label)

if ~isequal(size(X), size(Y))
    error('%s returned X and Y with different sizes: X=%s, Y=%s.', ...
        label, mat2str(size(X)), mat2str(size(Y)));
end

if ~isequal(size(X), [Nt, Ns]) && isequal(size(X), [Ns, Nt])
    warning('MOLECompare:TransposeGrid', '%s returned Ns x Nt. Transposing to Nt x Ns.', label);
    X = X.';
    Y = Y.';
end

if ~isequal(size(X), [Nt, Ns])
    error('Expected %s output size Nt x Ns = %d x %d, got %s.', ...
        label, Nt, Ns, mat2str(size(X)));
end

end

%% ========================================================================
function gridSizes = normalizeGridSizes(gridSizes)

if isempty(gridSizes)
    error('gridSizes cannot be empty. Use one [Ns Nt] row per run.');
end

if isvector(gridSizes)
    gridSizes = gridSizes(:).';
end

if size(gridSizes, 2) ~= 2
    error(['gridSizes must be an N-by-2 matrix with one [Ns Nt] row per run.\n', ...
           'Examples: gridSizes = [41 41]; or gridSizes = [21 21; 41 41; 81 81].']);
end

gridSizes = round(gridSizes);

if any(gridSizes(:) < 3)
    error('All Ns and Nt values must be at least 3.');
end

end

%% ========================================================================
function tf = shouldMakePlots(Ns, Nt, plotGridSize, makePlots)

if ~makePlots || isempty(plotGridSize)
    tf = false;
    return;
end

if isscalar(plotGridSize)
    tf = (Ns == plotGridSize && Nt == plotGridSize);
else
    tf = any(plotGridSize(:,1) == Ns & plotGridSize(:,2) == Nt);
end

end

%% ========================================================================
function out = validatestringCompat(value, allowed, name)

if isstring(value)
    value = char(value);
end
value = char(value);

matches = strcmpi(value, allowed);
if ~any(matches)
    error('%s must be one of: %s.', name, strjoinCompat(allowed, ', '));
end

out = allowed{find(matches, 1, 'first')};

end

%% ========================================================================
function s = strjoinCompat(parts, delimiter)

if exist('strjoin', 'builtin') || exist('strjoin', 'file')
    s = strjoin(parts, delimiter);
    return;
end

s = parts{1};
for ii = 2:numel(parts)
    s = [s delimiter parts{ii}]; %#ok<AGROW>
end

end

%% ========================================================================
function stem = sanitizeFileStem(textIn)

stem = regexprep(char(textIn), '[^A-Za-z0-9_\-]+', '_');
stem = regexprep(stem, '_+', '_');
stem = regexprep(stem, '^_+|_+$', '');
if isempty(stem)
    stem = 'grid';
end

end

%% ========================================================================
function applyCLim(climVals)

try
    clim(climVals);
catch
    caxis(climVals);
end

end
