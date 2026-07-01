function swan_length_grid_vanilla()
% Vanilla length grid (pure Laplace) for the SWAN boundaries in your script.
% Solves ΔX=0 and ΔY=0 on the unit (s,t) grid with your boundary curves,
% then plots the resulting grid.

%% ---------------- User settings ----------------
Ns = 41; Nt = 41;      % logical grid size
I = Ns; J = Nt;
hs = 1/(I-1); ht = 1/(J-1);

%% ---------------- Boundaries (from your code) ----------------
xb_fun = @(s) s;                  yb_fun = @(s) 0*s;
xt_fun = @(s) s;                  yt_fun = @(s) 1 - 3*s + 3*s.^2;
xl_fun = @(s) 0*s;                yl_fun = @(s) s;
xr_fun = @(s) 1 + 2*s - 2*s.^2;   yr_fun = @(s) s;

sb = linspace(0,1,I); sl = linspace(0,1,J).';
xb = xb_fun(sb);  yb = yb_fun(sb);
xt = xt_fun(sb);  yt = yt_fun(sb);
xl = xl_fun(sl);  yl = yl_fun(sl);
xr = xr_fun(sl);  yr = yr_fun(sl);

% Pack boundary arrays as (J x I): rows ~ t, cols ~ s
BX = nan(J,I); BY = nan(J,I);
BX(1,:) = xb;  BY(1,:) = yb;         % bottom (t=0)
BX(J,:) = xt;  BY(J,:) = yt;         % top    (t=1)
BX(:,1) = xl;  BY(:,1) = yl;         % left   (s=0)
BX(:,I) = xr;  BY(:,I) = yr;         % right  (s=1)
% corner consistency (already consistent, but ensure numeric equality)
BX(1,1)=xb(1);   BY(1,1)=yb(1);
BX(1,I)=xb(end); BY(1,I)=yb(end);
BX(J,1)=xt(1);   BY(J,1)=yt(1);
BX(J,I)=xt(end); BY(J,I)=yt(end);

%% ---------------- Laplace operators (Dirichlet interior) ----------------
nI = I-2; nJ = J-2;
eI = ones(nI,1); eJ = ones(nJ,1);
Ts = spdiags([-eI 2*eI -eI], [-1 0 1], nI, nI);  % s-Laplacian
Tt = spdiags([-eJ 2*eJ -eJ], [-1 0 1], nJ, nJ);  % t-Laplacian

% Scaling (ΔX/Δs^2 + ΔX/Δt^2 = 0) ⇒ Sylvester: (wy*Tt)U + U(wx*Ts) = F
wx = 1/(hs^2);
wy = 1/(ht^2);
Aop = wy * full(Tt);
Bop = wx * full(Ts);

% RHS from boundaries for X and Y
bFxL = BX(2:J-1, 1);   bFyL = BY(2:J-1, 1);
bFxR = BX(2:J-1, I);   bFyR = BY(2:J-1, I);
bFxB = BX(1, 2:I-1);   bFyB = BY(1, 2:I-1);
bFxT = BX(J, 2:I-1);   bFyT = BY(J, 2:I-1);

Fx = zeros(nJ,nI);  Fy = zeros(nJ,nI);
Fx(:,1)   = Fx(:,1)   + wx * bFxL;
Fx(:,end) = Fx(:,end) + wx * bFxR;
Fx(1,:)   = Fx(1,:)   + wy * bFxB;
Fx(end,:) = Fx(end,:) + wy * bFxT;

Fy(:,1)   = Fy(:,1)   + wx * bFyL;
Fy(:,end) = Fy(:,end) + wx * bFyR;
Fy(1,:)   = Fy(1,:)   + wy * bFyB;
Fy(end,:) = Fy(end,:) + wy * bFyT;

%% ---------------- Solve Sylvester systems ----------------
Ux = sylvester(Aop, Bop, Fx);   % (J-2) x (I-2)
Uy = sylvester(Aop, Bop, Fy);

% Insert interior into full grids
X = BX; Y = BY;
X(2:J-1, 2:I-1) = Ux;
Y(2:J-1, 2:I-1) = Uy;

%% ---------------- Plot vanilla length grid ----------------
figure('Color','w'); hold on; axis equal; box on;
col1 = [0.15 0.4 0.85];  % s-lines
col2 = [0.85 0.2 0.2];   % t-lines

% boundaries
plot(X( [1 end],:).', Y( [1 end],:).', 'k-', 'LineWidth', 1.25); % bottom & top
plot(X(:, [1 end]),    Y(:, [1 end]),    'k-', 'LineWidth', 1.25); % left & right

% interior grid lines (s = const; t = const)
for j = 1:J
    plot(X(j,:), Y(j,:), '-', 'Color', col1, 'LineWidth', 0.8);
end
for i = 1:I
    plot(X(:,i), Y(:,i), '-', 'Color', col2, 'LineWidth', 0.8);
end

title(sprintf('Vanilla length grid (Laplace)  Ns=%d, Nt=%d', Ns, Nt));
xlabel('x'); ylabel('y'); hold off;

%% ---------------- Optional: quick Jacobian check ----------------
%{
xs = (X(:,3:end)-X(:,1:end-2))/(2*hs);
ys = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
xt = (X(3:end,:)-X(1:end-2,:))/(2*ht);
yt = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);
Jc = xs(2:end-1,:).*yt(:,2:end-1) - xt(:,2:end-1).*ys(2:end-1,:);
fprintf('Min det(J) over cells: %.3e\n', min(Jc(:)));
%}
end
