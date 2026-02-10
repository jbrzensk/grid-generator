%% test_hardc_ref_grid
% Script to test gridGenerator
clc; clear; close all;

% Add Grids folder to path (where makeGridHardC.m lives)
addpath(fullfile(pwd,'Grids'));

%% --- Parameters ---
Ns = 23;       % # columns (s-direction)  (I)
Nt = 40;       % # rows    (t-direction)  (J)
N_basis = 0 ;  % # bubble basis functions

%% --- Grid function handle ---
grid_fun = @makeGridHardC;  % should accept (Ns-1, Nt-1) and return X,Y on boundary grid

%% --- Quick sanity check on grid_fun output sizes ---
try
    [Xc, Yc] = grid_fun(Ns-1, Nt-1);
catch ME
    error('grid_fun failed when called as grid_fun(Ns-1, Nt-1). Error:\n%s', ME.message);
end

fprintf('grid_fun returned size(Xc)=%s, size(Yc)=%s\n', mat2str(size(Xc)), mat2str(size(Yc)));

% Many implementations return (Nt x Ns). We want boundary arrays consistent with Nt,Ns.
% If yours returns (Ns x Nt), transpose them here.
if ~isequal(size(Xc), [Nt, Ns]) && isequal(size(Xc), [Ns, Nt])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.'; Yc = Yc.';
end

if ~isequal(size(Xc), [Nt, Ns])
    warning(['Expected grid_fun to return size (Nt x Ns) = (%d x %d). ' ...
             'Got %s. gridGenerator may error unless its boundary assembly matches this.'], ...
             Nt, Ns, mat2str(size(Xc)));
end

%% --- Run gridGenerator ---
try
    [C_opt, E_opt, Xopt, Yopt, Jc, Aopt, Bopt] = gridGenerator(grid_fun, Ns, Nt, N_basis);
catch ME
    disp(ME.getReport('extended'));
    error('gridGenerator failed. See error report above.');
end

%% --- Print results ---
fprintf('\nRESULTS:\n');
fprintf('  C_opt   = %.8f\n', C_opt);
fprintf('  E_opt   = %.6e (generalized solver NH energy)\n', E_opt);
fprintf('  min(Jc) = %.3e\n', min(Jc(:)));

if ~isempty(Aopt)
    fprintf('  ||A||_2 = %.3e, ||B||_2 = %.3e\n', norm(Aopt), norm(Bopt));
end

%% --- Plot final grid (if you want an extra plot) ---
figure; hold on; axis equal; box on;
for j=1:Nt
    plot(Xopt(j,:), Yopt(j,:), '-k');
end
for i=1:Ns
    plot(Xopt(:,i), Yopt(:,i), '-k');
end
title(sprintf('Final grid (C*=%.4f, Nbasis=%d)', C_opt, N_basis));
xlabel('X'); ylabel('Y');

%% --- Plot Jacobian heatmap ---
figure;
imagesc(Jc); axis image; colorbar;
title('Jacobian det J at cell centers');
xlabel('i'); ylabel('j');
