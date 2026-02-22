%% test_hardC_softmax_gridGenerator
% Script to test the *softmax-reclustering* version of gridGenerator on Hard-C
clc; clear; close all;

% Robust path setup (so pwd doesn't matter)
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir,'Grids'));

%% --- Parameters ---
Ns = 41;        % # columns (s-direction)
Nt = 41;        % # rows    (t-direction)

% NOTE: In the modified gridGenerator, this is interpreted as:
%   N_basis = number of SGD steps for softmax reclustering
N_basis = 1000;  % try 100..1000

%% --- Grid function handle ---
grid_fun = @makeGridHardC;   % must return [X,Y] of size Nt x Ns

%% --- Sanity check (optional but useful) ---
[Xc, Yc] = grid_fun(Ns-1, Nt-1);
fprintf('grid_fun returned size(Xc)=%s, size(Yc)=%s\n', mat2str(size(Xc)), mat2str(size(Yc)));

if ~isequal(size(Xc), [Nt, Ns]) && isequal(size(Xc), [Ns, Nt])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.'; Yc = Yc.';
end
if ~isequal(size(Xc), [Nt, Ns])
    warning('Expected (Nt x Ns) = (%d x %d), got %s. Proceeding anyway...', Nt, Ns, mat2str(size(Xc)));
end

%% --- Run the optimization (plots are produced by gridGenerator itself) ---
% The modified gridGenerator returns and plots:
%   - Physical grid
%   - Jacobian heatmap
%   - Reference (u,v) grid
%   - Overlay reference vs physical
[C_opt, E_opt, Xopt, Yopt, Jc, pU_opt, pV_opt, Uopt, Vopt] = gridGenerator(grid_fun, Ns, Nt, N_basis);

%% --- Print results ---
fprintf('\nRESULTS:\n');
fprintf('  C_opt   = %.8f\n', C_opt);
fprintf('  E_opt   = %.6e\n', E_opt);
fprintf('  min(Jc) = %.3e\n', min(Jc(:)));

%% --- Extra plots (optional): spacing distributions in the reclustered reference grid ---
% pU_opt and pV_opt are logits; convert to spacings
dU = softmax_local(pU_opt(:));
dV = softmax_local(pV_opt(:));
dU = dU / sum(dU);
dV = dV / sum(dV) * C_opt;

figure; subplot(2,1,1);
plot(dU,'-o'); grid on;
title('Reference spacings dU (sum=1)'); xlabel('i'); ylabel('\Delta U_i');

subplot(2,1,2);
plot(dV,'-o'); grid on;
title(sprintf('Reference spacings dV (sum=C=%.4f)', C_opt)); xlabel('j'); ylabel('\Delta V_j');

%% --- Local softmax (uses built-in if available) ---
function y = softmax_local(x)
x = x(:);
if exist('softmax','file') == 2
    y = softmax(x);
else
    x = x - max(x);
    ex = exp(x);
    y = ex / sum(ex);
end
end
