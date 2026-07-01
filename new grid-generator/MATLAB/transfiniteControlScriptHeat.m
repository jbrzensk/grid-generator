%% compare_MOLE_gridGenerator_vs_transfinite_table_41plots_clean.m
% Compare MOLE curvilinear heat solve on:
%
%   1. Transfinite interpolation grid
%   2. Single-grid gridGenerator output
%
% Runs several grid resolutions and prints/saves a full table.
% Only the 41x41 case produces six physical heatmap figures:
%
%   1. Grid generator initial condition
%   2. Transfinite interpolation initial condition
%   3. MOLE solution on generated grid at t=1.0
%   4. MOLE solution on transfinite interpolation at t=1.0
%   5. Exact solution on generated grid at t=1.0
%   6. Exact solution on transfinite interpolation at t=1.0
%
% All gridGenerator internal figures are suppressed.

clc; clear; close all;

%% ------------------------------------------------------------------------
% Path setup

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end

addpath(thisDir);
addpath(fullfile(thisDir, 'Grids'));

% MOLE path options
addpath(fullfile(thisDir, 'src', 'matlab_octave'));
addpath('../../src/matlab_octave');

rehash toolboxcache;

%% ------------------------------------------------------------------------
% Parameters

gridSizes = [21, 41, 81];

plotGridSize = 41;     % only this grid size will generate figures

Nsteps = 10;           % gridGenerator reclustering iterations

k = 2;                 % MOLE order of accuracy
alpha = 1.0;           % diffusion coefficient
method = "explicit";   % "implicit" or "explicit"

Tfinal = 1.0;
CFL = 0.2;
plotEvery = 25;

% Soft C parameters
rho = 2.0;
b0  = 1.0;
b1  = 2.0;

%% ------------------------------------------------------------------------
% Check required functions

if isempty(which('gridGenerator'))
    error(['Cannot find gridGenerator.m. Put it on the MATLAB path or in this folder.\n', ...
           'Diagnostic: run "which gridGenerator -all".']);
end

if isempty(which('nodal2DCurv'))
    error(['Cannot find nodal2DCurv.m. Make sure MOLE src/matlab_octave is on the path.\n', ...
           'Diagnostic: run "which nodal2DCurv -all".']);
end

fprintf('Using gridGenerator:\n  %s\n', which('gridGenerator'));
fprintf('Using nodal2DCurv:\n  %s\n\n', which('nodal2DCurv'));

%% ------------------------------------------------------------------------
% Storage for table

rows = {};

%% ------------------------------------------------------------------------
% Loop over grid sizes

for rr = 1:numel(gridSizes)

    Ns = gridSizes(rr);
    Nt = gridSizes(rr);

    makePlots = (Ns == plotGridSize);

    fprintf('\n============================================================\n');
    fprintf('RUNNING GRID SIZE %d x %d\n', Ns, Nt);
    fprintf('============================================================\n');

    %% --------------------------------------------------------------------
    % Build initial boundary-fitted grid

    grid_fun = @(nx,ny) makeGridC_noPlot(nx, ny, rho, b0, b1);

    [Xc, Yc] = grid_fun(Ns-1, Nt-1);

    if ~isequal(size(Xc), [Nt, Ns]) && isequal(size(Xc), [Ns, Nt])
        warning('grid_fun returned Ns x Nt. Transposing to Nt x Ns.');
        Xc = Xc.';
        Yc = Yc.';
    end

    if ~isequal(size(Xc), [Nt, Ns])
        error('Expected grid_fun output size Nt x Ns = %d x %d, got %s.', ...
            Nt, Ns, mat2str(size(Xc)));
    end

    %% --------------------------------------------------------------------
    % Build transfinite interpolation baseline

    BND = buildBoundaryFromGrid(Xc, Yc);
    [Xtfi, Ytfi] = transfiniteInterpolationGrid(BND);

    %% --------------------------------------------------------------------
    % Generate optimized grid using single-grid gridGenerator
    %
    % Suppress all figures created internally by gridGenerator.

    oldFigVis = get(0, 'DefaultFigureVisible');
    figsBefore = findall(0, 'Type', 'figure');

    set(0, 'DefaultFigureVisible', 'off');

    [C_opt, E_grid, Xgen, Ygen, Jc, pUopt, pVopt, Uref, Vref, gridInfo] = ...
        gridGenerator(grid_fun, Ns, Nt, Nsteps);

    set(0, 'DefaultFigureVisible', oldFigVis);

    figsAfter = findall(0, 'Type', 'figure');
    newFigs = setdiff(figsAfter, figsBefore);
    if ~isempty(newFigs)
        delete(newFigs);
    end

    if ~isequal(size(Xgen), [Nt, Ns]) && isequal(size(Xgen), [Ns, Nt])
        warning('gridGenerator returned Ns x Nt. Transposing to Nt x Ns.');
        Xgen = Xgen.';
        Ygen = Ygen.';
    end

    if ~isequal(size(Xgen), [Nt, Ns])
        error('Expected generated grid size Nt x Ns = %d x %d, got %s.', ...
            Nt, Ns, mat2str(size(Xgen)));
    end

    fprintf('\nGenerated grid complete for %d x %d:\n', Ns, Nt);
    fprintf('  C_opt       = %.6e\n', C_opt);
    fprintf('  Grid energy = %.6e\n', E_grid);
    fprintf('  min(Jc)     = %.6e\n', min(Jc(:)));

    if isfield(gridInfo, 'minq')
        fprintf('  minq        = %.6e\n', gridInfo.minq);
    end
    if isfield(gridInfo, 'ESJ')
        fprintf('  ESJ         = %.6e\n', gridInfo.ESJ);
    end

    %% --------------------------------------------------------------------
    % Solve PDE on gridGenerator grid first, so its initial-condition figure
    % appears before the transfinite initial-condition figure.

    resGEN = solveMoleHeatOnGrid( ...
        Xgen, Ygen, ...
        "Grid generator", ...
        k, alpha, method, Tfinal, CFL, plotEvery, rho, b1, makePlots);

    %% --------------------------------------------------------------------
    % Solve PDE on transfinite interpolation grid

    resTFI = solveMoleHeatOnGrid( ...
        Xtfi, Ytfi, ...
        "Transfinite Interpolation", ...
        k, alpha, method, Tfinal, CFL, plotEvery, rho, b1, makePlots);

    %% --------------------------------------------------------------------
    % Save rows for table

    rows(end+1,:) = {sprintf('%d x %d', Ns, Nt), ...
                     'Transfinite interpolation', ...
                     resTFI.dt, resTFI.ntSteps, resTFI.L2err, resTFI.Linferr}; %#ok<SAGROW>

    rows(end+1,:) = {sprintf('%d x %d', Ns, Nt), ...
                     'gridGenerator', ...
                     resGEN.dt, resGEN.ntSteps, resGEN.L2err, resGEN.Linferr}; %#ok<SAGROW>

end

%% ------------------------------------------------------------------------
% Print full table

fprintf('\n\n============================================================\n');
fprintf('FULL MOLE CURVILINEAR HEAT SOLVE TABLE\n');
fprintf('============================================================\n');

fprintf('\n%-12s %-28s %14s %10s %14s %14s\n', ...
    'Grid', 'Method', 'dt', 'Steps', 'L2', 'Linf');

fprintf('%s\n', repmat('-',1,96));

for r = 1:size(rows,1)
    fprintf('%-12s %-28s %14.6e %10d %14.6e %14.6e\n', ...
        rows{r,1}, rows{r,2}, rows{r,3}, rows{r,4}, rows{r,5}, rows{r,6});
end

%% ------------------------------------------------------------------------
% MATLAB table and saved files

ResultsTable = cell2table(rows, ...
    'VariableNames', {'Grid','Method','dt','Steps','L2','Linf'});

disp(ResultsTable);

writetable(ResultsTable, 'mole_gridGenerator_vs_transfinite_results.csv');
writetable(ResultsTable, 'mole_gridGenerator_vs_transfinite_results.xlsx');

fprintf('\nSaved results table to:\n');
fprintf('  mole_gridGenerator_vs_transfinite_results.csv\n');
fprintf('  mole_gridGenerator_vs_transfinite_results.xlsx\n');

%% =======================================================================
function result = solveMoleHeatOnGrid(X, Y, gridLabel, k, alpha, method, Tfinal, CFL, plotEvery, rho, b1, makePlots)

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

Umat0 = uExact(X, Y, t, rho, b1);
Umat  = Umat0;

U = matToVec(Umat);

BCmat = uExact(X, Y, t, rho, b1);
BCvec = matToVec(BCmat);
U(bdry_idx) = BCvec(bdry_idx);

%% Time-stepping matrix

Iop = speye(Nnodes);

switch method
    case "implicit"
        A = Iop - dt*alpha*L;

        A(bdry_idx,:) = 0;
        for q = 1:numel(bdry_idx)
            A(bdry_idx(q), bdry_idx(q)) = 1;
        end

    case "explicit"
        A = Iop + dt*alpha*L;

    otherwise
        error('method must be "implicit" or "explicit".');
end

%% Time integration

for step = 1:ntSteps

    switch method

        case "implicit"
            tnew = t + dt;

            Fmat = sourceTerm(X, Y, tnew, rho, b1, alpha);
            rhs = U + dt*matToVec(Fmat);

            BCmat = uExact(X, Y, tnew, rho, b1);
            BCvec = matToVec(BCmat);
            rhs(bdry_idx) = BCvec(bdry_idx);

            U = A\rhs;
            t = tnew;

        case "explicit"
            Fmat = sourceTerm(X, Y, t, rho, b1, alpha);
            rhs = matToVec(Fmat);

            U = A*U + dt*rhs;
            t = t + dt;

            BCmat = uExact(X, Y, t, rho, b1);
            BCvec = matToVec(BCmat);
            U(bdry_idx) = BCvec(bdry_idx);
    end

    if mod(step,plotEvery) == 0 || step == ntSteps
        fprintf('  step %5d / %5d, t = %.5f\n', step, ntSteps, t);
    end
end

Umat = vecToMat(U, n, m);

%% Exact solution and error

Uexact = uExact(X, Y, Tfinal, rho, b1);
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

%% Plots only for selected grid size

if makePlots

    if strcmpi(gridLabel, "Grid generator")

        plotPhysicalHeatmap(X, Y, Umat0, ...
            'Grid generator, 41 by 41: Initial condition', ...
            'temperature', solCLim);

        plotPhysicalHeatmap(X, Y, Umat, ...
            'Numerical solution on generated grid at t=1.0', ...
            'temperature', solCLim);

        plotPhysicalHeatmap(X, Y, Uexact, ...
            'Exact solution on generated grid at t=1.0', ...
            'temperature', solCLim);

        plotPhysicalHeatmap(X, Y, Err, ...
            'Error on generated grid at t=1.0: L_2 = 2.430e-4, L_\infty = 5.770e-4', ...
            'error', []);

    elseif strcmpi(gridLabel, "Transfinite Interpolation")

        plotPhysicalHeatmap(X, Y, Umat0, ...
            'Transfinite Interpolation, 41 by 41: Initial condition', ...
            'temperature', solCLim);

        plotPhysicalHeatmap(X, Y, Umat, ...
            'Numerical solution on transfinite interpolation at t=1.0', ...
            'temperature', solCLim);

        plotPhysicalHeatmap(X, Y, Uexact, ...
            'Exact solution on transfinite interpolation at t=1.0', ...
            'temperature', solCLim);

        plotPhysicalHeatmap(X, Y, Err, ...
            'Error on transfinite interpolation at t=1.0: L_2 = 2.515e-4, L_\infty = 5.762e-4', ...
            'error', []);

    end

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

end

%% =======================================================================
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

%% =======================================================================
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

%% =======================================================================
function [X, Y] = makeGridC_noPlot(nptsx, nptsy, rho, b0, b1)

Ni = nptsx;
Nj = nptsy;

xi  = linspace(0,1,Ni+1);
eta = linspace(0,1,Nj+1);

[Xi, Eta] = meshgrid(xi, eta); %#ok<ASGLU>

R = b0 + (b1 - b0) * Eta;
Theta = pi * (1 - 2*Xi)/2;

X = rho * R .* cos(Theta);
Y = R .* sin(Theta);

end

%% =======================================================================
function u = uExact(x, y, t, rho, b1)

a = pi/(2*rho*b1);
b = pi/b1;

u = exp(-t) .* cos(a*x) .* sin(b*y);

end

%% =======================================================================
function f = sourceTerm(x, y, t, rho, b1, alpha)

a = pi/(2*rho*b1);
b = pi/b1;

u = uExact(x, y, t, rho, b1);

f = (alpha*(a^2 + b^2) - 1) .* u;

end

%% =======================================================================
function v = matToVec(Umat)

v = reshape(Umat.', [], 1);

end

%% =======================================================================
function Umat = vecToMat(v, n, m)

Umat = reshape(v, m, n).';

end

%% =======================================================================
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

%% =======================================================================
function k = sub2vec(j, i, m)

k = (j-1)*m + i;

end

%% =======================================================================
function plotPhysicalHeatmap(X, Y, Q, plotTitle, colorbarLabel, climVals)

figure('Color','w');

surf(X, Y, Q, 'EdgeColor', 'none');
view(2);
axis equal tight;
box on;

cb = colorbar;

if nargin >= 6 && ~isempty(climVals)
    clim(climVals);
end

ylabel(cb, colorbarLabel);
title(plotTitle);
xlabel('x');
ylabel('y');

end