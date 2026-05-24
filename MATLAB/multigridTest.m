%% test_compactAC_multigrid
clc; clear; close all;

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir,'Grids'));

%% --- Parameters ---
Ns = 81;
Nt = 81;
N_outer = 0000;   % number of SGD/equalization steps

%% --- Grid function handle ---
grid_fun = @makeGridA;  % change to @makeGridC if desired

%% --- Sanity check on grid_fun output sizes ---
[Xc, Yc] = grid_fun(Ns-1, Nt-1);
fprintf('grid_fun returned size(Xc)=%s, size(Yc)=%s\n', ...
    mat2str(size(Xc)), mat2str(size(Yc)));

if ~isequal(size(Xc), [Nt, Ns]) && isequal(size(Xc), [Ns, Nt])
    warning('grid_fun returned (Ns x Nt). Transposing to (Nt x Ns).');
    Xc = Xc.';
    Yc = Yc.';
end

if ~isequal(size(Xc), [Nt, Ns])
    warning('Expected (Nt x Ns)=(%d x %d), got %s. Proceeding anyway...', ...
        Nt, Ns, mat2str(size(Xc)));
end

%% --- Run compact [0,Inf]^2 A,C multigrid generator ---
[A_opt, C_opt, E_opt, Xopt, Yopt, Jc, ...
    pUopt, pVopt, pBopt, pTopt, pLopt, pRopt, ...
    Uopt, Vopt, info] = gridGeneratorMultigrid( ...
        grid_fun, Ns, Nt, N_outer);

%% --- Print results ---
fprintf('\nRESULTS:\n');
fprintf('  A_opt   = %s\n', compactNumToStr(A_opt));
fprintf('  C_opt   = %s\n', compactNumToStr(C_opt));
fprintf('  E_opt   = %.6e\n', E_opt);
fprintf('  min(Jc) = %.3e\n', min(Jc(:)));
fprintf('  max(Jc) = %.3e\n', max(Jc(:)));

if isfield(info,'accepted')
    fprintf('  accepted steps = %d\n', info.accepted);
end

if isfield(info,'thetaW')
    fprintf('  thetaW  = %.6f\n', info.thetaW);
end

if isfield(info,'minq')
    fprintf('  minq    = %.6e\n', info.minq);
end

if isfield(info,'ESJ')
    fprintf('  ESJ     = %.6e\n', info.ESJ);
end

%% --- Reference spacings from pUopt,pVopt ---
refHeight = 1.0;
if isfield(info,'refHeight')
    refHeight = info.refHeight;
end

dU = expStable(pUopt(:));
dU = dU / sum(dU);

dV = expStable(pVopt(:));
dV = dV / sum(dV) * refHeight;

figure;
subplot(2,1,1);
plot(dU,'-o'); grid on;
title('Reference spacings dU, sum = 1');
xlabel('i'); ylabel('\Delta U_i');

subplot(2,1,2);
plot(dV,'-o'); grid on;
title(sprintf('Reference spacings dV, sum = refHeight = %.4f', refHeight));
xlabel('j'); ylabel('\Delta V_j');

%% --- Physical boundary spacings ---
dB = expStable(pBopt(:)); dB = dB / sum(dB);
dT = expStable(pTopt(:)); dT = dT / sum(dT);
dL = expStable(pLopt(:)); dL = dL / sum(dL);
dR = expStable(pRopt(:)); dR = dR / sum(dR);

figure;
subplot(2,2,1);
plot(dB,'-o'); grid on;
title('Bottom boundary spacings');
xlabel('i'); ylabel('\Delta s_B');

subplot(2,2,2);
plot(dT,'-o'); grid on;
title('Top boundary spacings');
xlabel('i'); ylabel('\Delta s_T');

subplot(2,2,3);
plot(dL,'-o'); grid on;
title('Left boundary spacings');
xlabel('j'); ylabel('\Delta s_L');

subplot(2,2,4);
plot(dR,'-o'); grid on;
title('Right boundary spacings');
xlabel('j'); ylabel('\Delta s_R');

%% --- Plot physical grid ---
figure; hold on; axis equal; box on;
for jj = 1:Nt
    plot(Xopt(jj,:), Yopt(jj,:), '-k');
end
for ii = 1:Ns
    plot(Xopt(:,ii), Yopt(:,ii), '-k');
end
title(sprintf('Physical grid, A=%s, C=%s', ...
    compactNumToStr(A_opt), compactNumToStr(C_opt)));
xlabel('X'); ylabel('Y');

%% --- Plot reference grid ---
figure; hold on; axis equal; box on;
for jj = 1:Nt
    plot(Uopt(jj,:), Vopt(jj,:), '-b');
end
for ii = 1:Ns
    plot(Uopt(:,ii), Vopt(:,ii), '-b');
end
title('Reference grid after simultaneous equalization');
xlabel('U'); ylabel('V');

%% --- Plot Jacobian ---
figure;
imagesc(Jc);
axis image;
colorbar;
title('Jacobian Jc');
xlabel('i'); ylabel('j');

%% --- History plots ---
if isfield(info,'Ehist') && ~isempty(info.Ehist)
    figure;
    plot(info.Ehist,'-o'); grid on;
    title('Energy history');
    xlabel('iteration'); ylabel('E');
end

if isfield(info,'minJhist') && ~isempty(info.minJhist)
    figure;
    plot(info.minJhist,'-o'); grid on;
    title('Minimum Jacobian history');
    xlabel('iteration'); ylabel('min J');
end

if isfield(info,'negHist') && ~isempty(info.negHist)
    figure;
    plot(info.negHist,'-o'); grid on;
    title('Negative-cell count history');
    xlabel('iteration'); ylabel('# cells with J <= 0');
end

if isfield(info,'vlogJhist') && ~isempty(info.vlogJhist)
    figure;
    plot(info.vlogJhist,'-o'); grid on;
    title('Var(log J_+) history');
    xlabel('iteration'); ylabel('Var(log J_+)');
end

if isfield(info,'stepHist') && ~isempty(info.stepHist)
    figure;
    plot(info.stepHist,'-o'); grid on;
    title('Accepted step history');
    xlabel('iteration'); ylabel('step size');
end

if isfield(info,'ghist') && ~isempty(info.ghist)
    figure;
    plot(info.ghist,'-o'); grid on;
    title('Sliver metric history');
    xlabel('iteration'); ylabel('g');
end

%% --- Optional: Stage 1 A,C diagnostics ---
if isfield(info,'AC')
    disp('Stage 1 compact A,C search info:');
    disp(info.AC);
end

%% ========================================================================
%% Local helpers

function y = expStable(x)
    x = x(:);
    x = x - max(x);
    y = exp(x);
end

function s = compactNumToStr(x)
    if isinf(x)
        if x > 0
            s = 'Inf';
        else
            s = '-Inf';
        end
    elseif x == 0
        s = '0';
    elseif isnan(x)
        s = 'NaN';
    else
        s = sprintf('%.6g', x);
    end
end