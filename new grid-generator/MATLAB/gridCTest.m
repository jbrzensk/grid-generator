%% test_softC
clc; clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir,'Grids'));

%% --- Parameters ---
Ns = 41;
Nt = 21;
N_outer = 0010;   % number of ACCEPTED SGD steps (with filter)

%% --- Grid function handle ---
grid_fun = @makeGridChevron;  % <-- change to @makeGridC if you really want the C-grid

%% --- Sanity check on grid_fun output sizes ---
[Xc, Yc] = grid_fun(Ns-1, Nt-1);
fprintf('grid_fun returned size(Xc)=%s, size(Yc)=%s\n', mat2str(size(Xc)), mat2str(size(Yc)));

if ~isequal(size(Xc), [Nt, Ns]) && isequal(size(Xc), [Ns, Nt])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.'; Yc = Yc.';
end
if ~isequal(size(Xc), [Nt, Ns])
    warning('Expected (Nt x Ns)=(%d x %d), got %s. Proceeding anyway...', Nt, Ns, mat2str(size(Xc)));
end

%% --- Run ---
[C_opt, E_opt, Xopt, Yopt, Jc, pUopt, pVopt, Uopt, Vopt, info] = gridGenerator(grid_fun, Ns, Nt, N_outer);

%% --- Print results ---
fprintf('\nRESULTS:\n');
fprintf('  C_opt   = %.6f\n', C_opt);
fprintf('  E_opt   = %.6e\n', E_opt);
fprintf('  min(Jc) = %.3e\n', min(Jc(:)));

if isfield(info,'accepted')
    fprintf('  accepted steps = %d\n', info.accepted);
end

%% --- Plot reference spacings from pUopt,pVopt ---
dU = exp(pUopt(:)); dU = dU / sum(dU);               % sum=1
dV = exp(pVopt(:)); dV = dV / sum(dV) * C_opt;       % sum=C_opt

figure;
subplot(2,1,1);
plot(dU,'-o'); grid on;
title('Reference spacings dU (sum=1)'); xlabel('i'); ylabel('\Delta U_i');

subplot(2,1,2);
plot(dV,'-o'); grid on;
title(sprintf('Reference spacings dV (sum=C=%.4f)', C_opt)); xlabel('j'); ylabel('\Delta V_j');

%% --- History plots (if present) ---
if isfield(info,'Ehist') && ~isempty(info.Ehist)
    figure; plot(info.Ehist,'-o'); grid on;
    title('Energy history'); xlabel('accepted step'); ylabel('E');
end
if isfield(info,'ghist') && ~isempty(info.ghist)
    figure; plot(info.ghist,'-o'); grid on;
    title('Sliver metric history'); xlabel('accepted step'); ylabel('g');
end