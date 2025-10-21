%% swan_optimize_B_hyperelastic.m
% Find the best B (with A normalized to 1) that minimizes the hyperelastic
% energy (Eq. 4) on the Swan grid built via the reference-length functional (Eq. 12).
% The grid for a given B is computed fast via a Sylvester solve (separable case).

clear; clc;

%% --- User settings --------------------------------------------------------
% Logical mesh sizes (choose roughly I:J ~ 1:B for nicely shaped cells)
I = 81;                       % nodes along s (columns)
J = 81;                       % nodes along t (rows)
A = 1.0;                       % normalize A to 1
B_bracket = [0.10, 2.00];      % search range for B (adjust as desired)

% Lamé parameters for Eq.(4) (scale doesn't matter; relative weighting does)
mu = 1.0;
lambda = 1.0;

% Boundary curves for the "Swan" (CCW, s in [0,1])
xb_fun = @(s) s;                  yb_fun = @(s) 0*s;                 % bottom
xt_fun = @(s) s;                  yt_fun = @(s) 1 - 3*s + 3*s.^2;    % top
xl_fun = @(s) 0*s;                yl_fun = @(s) s;                   % left
xr_fun = @(s) 1 + 2*s - 2*s.^2;   yr_fun = @(s) s;                   % right

% Build boundary arrays on the logical grid
sb = linspace(0,1,I); sl = linspace(0,1,J).';
xb = xb_fun(sb);  yb = yb_fun(sb);
xt = xt_fun(sb);  yt = yt_fun(sb);
xl = xl_fun(sl);  yl = yl_fun(sl);
xr = xr_fun(sl);  yr = yr_fun(sl);

% Pack once for reuse
BND.X = nan(J,I); BND.Y = nan(J,I);
BND.X(1,:) = xb;  BND.Y(1,:) = yb;
BND.X(J,:) = xt;  BND.Y(J,:) = yt;
BND.X(:,1) = xl;  BND.Y(:,1) = yl;
BND.X(:,I) = xr;  BND.Y(:,I) = yr;
BND.X(1,1)=xb(1);   BND.Y(1,1)=yb(1);
BND.X(1,I)=xb(end); BND.Y(1,I)=yb(end);
BND.X(J,1)=xt(1);   BND.Y(J,1)=yt(1);
BND.X(J,I)=xt(end); BND.Y(J,I)=yt(end);

% Logical spacings on [0,1]x[0,1]
hs = 1/(I-1); ht = 1/(J-1);

%% --- Objective: total hyperelastic energy for a given B -------------------
obj = @(B) total_energy_for_B(B, A, BND, I, J, hs, ht, mu, lambda);

% Quick check at a few B values (optional verbose probe)
% for Btest = linspace(B_bracket(1),B_bracket(2),5)
%     [E, minJ] = obj(Btest);
%     fprintf('B=%.3f  E=%.6e  minJ=%.3e\n', Btest, E, minJ);
% end

%% --- 1D minimize over B ---------------------------------------------------
opts = optimset('TolX',1e-3,'Display','iter');
[B_opt, E_opt] = fminbnd(@(B) obj(B), B_bracket(1), B_bracket(2), opts);

% Build final grid for best B
[~, ~, Xopt, Yopt, Jc] = obj(B_opt);

fprintf('\nBest B ≈ %.6f\nTotal energy ≈ %.6e\nMin Jacobian ≈ %.3e\n', B_opt, E_opt, min(Jc(:)));

%% --- Plot grid and Jacobian ----------------------------------------------
figure; hold on; axis equal; box on;
for j=1:J, plot(Xopt(j,:),Yopt(j,:),'-'); end
for i=1:I, plot(Xopt(:,i),Yopt(:,i),'-'); end
title(sprintf('Swan grid at B*=%.3f (A=1).', B_opt)); xlabel('x'); ylabel('y');

figure; imagesc(Jc); axis image; colorbar;
title('Jacobian det J at cell centers'); xlabel('i'); ylabel('j');

%% ======================= Helper functions ================================

function [E, minJ, X, Y, Jc] = total_energy_for_B(B, A, BND, I, J, hs, ht, mu, lambda)
    % Build grid for (A,B) via separable Eq.(12) and compute hyperelastic energy (Eq.4)

    % Fast Sylvester solve (constant A,B): weights
    wx = 1/(A^2*hs^2);                      % s-direction
    wy = 1/(B^2*ht^2);                      % t-direction

    % 1D Dirichlet second-difference (interior)
    nI = I-2; nJ = J-2;
    eI = ones(nI,1); eJ = ones(nJ,1);
    Ts = spdiags([-eI 2*eI -eI], [-1 0 1], nI, nI);
    Tt = spdiags([-eJ 2*eJ -eJ], [-1 0 1], nJ, nJ);
    Aop = wy * full(Tt);
    Bop = wx * full(Ts);

    % RHS from boundaries
    X = BND.X; Y = BND.Y;
    Fx = zeros(nJ,nI);  Fy = zeros(nJ,nI);
    Fx(:,1)   = Fx(:,1)   + wx * X(2:J-1, 1);   Fy(:,1)   = Fy(:,1)   + wx * Y(2:J-1, 1);
    Fx(:,end) = Fx(:,end) + wx * X(2:J-1, I);   Fy(:,end) = Fy(:,end) + wx * Y(2:J-1, I);
    Fx(1,:)   = Fx(1,:)   + wy * X(1, 2:I-1);   Fy(1,:)   = Fy(1,:)   + wy * Y(1, 2:I-1);
    Fx(end,:) = Fx(end,:) + wy * X(J, 2:I-1);   Fy(end,:) = Fy(end,:) + wy * Y(J, 2:I-1);

    % Solve for interior
    Ux = sylvester(Aop, Bop, Fx);
    Uy = sylvester(Aop, Bop, Fy);
    X(2:J-1,2:I-1) = Ux;
    Y(2:J-1,2:I-1) = Uy;

    % Metrics at cell centers (central diffs on unit logical grid)
    x_s = (X(:,3:end)-X(:,1:end-2))/(2*hs);
    y_s = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
    x_t = (X(3:end,:)-X(1:end-2,:))/(2*ht);
    y_t = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);

    xs = x_s(2:end-1,:);  ys = y_s(2:end-1,:);
    xt = x_t(:,2:end-1);  yt = y_t(:,2:end-1);

    % Jacobian and Right Cauchy-Green tensor pieces
    Jc = xs.*yt - xt.*ys;                   % det(F)
    minJ = min(Jc(:));

    % Hard barrier if any folding
    if minJ <= 0
        E = 1e12 + 1e8 * sum(Jc(:) <= 0);
        return
    end

    % 2D hyperelastic energy density (Eq.4 adapted with d=2)
    % C = F^T F, so tr(C) = xs.^2 + ys.^2 + xt.^2 + yt.^2
    trC = xs.^2 + ys.^2 + xt.^2 + yt.^2;
    logJ = log(Jc);
    d = 2;  % spatial dimension
    W = 0.5*mu*(trC - d) - mu*logJ + 0.5*lambda*(logJ.^2);

    % Integrate over logical domain (sum over cell centers)
    cellArea = hs*ht;
    E = cellArea * sum(W(:));
end
