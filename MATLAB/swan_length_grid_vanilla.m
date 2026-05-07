function swan_length_grid_vanilla()
% Vanilla length grid for a horseshoe / C-grid.
% Solves Laplace equations:
%    ΔX = 0,   ΔY = 0
% on the unit logical domain (s,t) in [0,1]x[0,1],
% with Dirichlet boundary data taken from the C-grid mapping.

%% ---------------- User settings ----------------
Ns = 21; 
Nt = 21;
I = Ns; 
J = Nt;

hs = 1/(I-1); 
ht = 1/(J-1);

%% ---------------- Horseshoe / C-grid boundaries ----------------
rho = 2.0;
b0  = 1.0;   % inner radius
b1  = 2.0;   % outer radius

% Boundary parameter values
sb = linspace(0,1,I);      % along bottom/top
sl = linspace(0,1,J).';    % along left/right (column vector)

% Boundary functions
xb_fun = @(s) rho*b0*cos(pi*(1 - 2*s)/2);
yb_fun = @(s)     b0*sin(pi*(1 - 2*s)/2);

xt_fun = @(s) rho*b1*cos(pi*(1 - 2*s)/2);
yt_fun = @(s)     b1*sin(pi*(1 - 2*s)/2);

xl_fun = @(t) 0*t;
yl_fun = @(t) b0 + (b1-b0)*t;

xr_fun = @(t) 0*t;
yr_fun = @(t) -(b0 + (b1-b0)*t);

% Evaluate boundaries
xb = xb_fun(sb);   yb = yb_fun(sb);
xt = xt_fun(sb);   yt = yt_fun(sb);
xl = xl_fun(sl);   yl = yl_fun(sl);
xr = xr_fun(sl);   yr = yr_fun(sl);

%% ---------------- Pack boundary arrays as (J x I) ----------------
% rows ~ t, cols ~ s
BX = nan(J,I); 
BY = nan(J,I);

BX(1,:) = xb;   BY(1,:) = yb;   % bottom (t=0)  inner arc
BX(J,:) = xt;   BY(J,:) = yt;   % top    (t=1)  outer arc
BX(:,1) = xl;   BY(:,1) = yl;   % left   (s=0)  upper radial side
BX(:,I) = xr;   BY(:,I) = yr;   % right  (s=1)  lower radial side

% Enforce exact corner consistency
BX(1,1) = xb(1);     BY(1,1) = yb(1);
BX(1,I) = xb(end);   BY(1,I) = yb(end);
BX(J,1) = xt(1);     BY(J,1) = yt(1);
BX(J,I) = xt(end);   BY(J,I) = yt(end);

%% ---------------- Laplace operators ----------------
nI = I-2; 
nJ = J-2;

eI = ones(nI,1); 
eJ = ones(nJ,1);

Ts = spdiags([-eI 2*eI -eI], [-1 0 1], nI, nI);  % s-direction
Tt = spdiags([-eJ 2*eJ -eJ], [-1 0 1], nJ, nJ);  % t-direction

wx = 1/(hs^2);
wy = 1/(ht^2);

Aop = wy * full(Tt);
Bop = wx * full(Ts);

%% ---------------- RHS from boundaries ----------------
bFxL = BX(2:J-1, 1);   bFyL = BY(2:J-1, 1);
bFxR = BX(2:J-1, I);   bFyR = BY(2:J-1, I);
bFxB = BX(1, 2:I-1);   bFyB = BY(1, 2:I-1);
bFxT = BX(J, 2:I-1);   bFyT = BY(J, 2:I-1);

Fx = zeros(nJ,nI);  
Fy = zeros(nJ,nI);

Fx(:,1)   = Fx(:,1)   + wx * bFxL;
Fx(:,end) = Fx(:,end) + wx * bFxR;
Fx(1,:)   = Fx(1,:)   + wy * bFxB;
Fx(end,:) = Fx(end,:) + wy * bFxT;

Fy(:,1)   = Fy(:,1)   + wx * bFyL;
Fy(:,end) = Fy(:,end) + wx * bFyR;
Fy(1,:)   = Fy(1,:)   + wy * bFyB;
Fy(end,:) = Fy(end,:) + wy * bFyT;

%% ---------------- Solve Sylvester systems ----------------
Ux = sylvester(Aop, Bop, Fx);
Uy = sylvester(Aop, Bop, Fy);

%% ---------------- Insert interior ----------------
X = BX; 
Y = BY;

X(2:J-1, 2:I-1) = Ux;
Y(2:J-1, 2:I-1) = Uy;

%% ---------------- Plot grid ----------------
figure('Color','w'); 
hold on; 
axis equal; 
box on;

col1 = [0.15 0.4 0.85];   % row lines
col2 = [0.85 0.2 0.2];    % column lines

% boundary curves
plot(X([1 end],:).', Y([1 end],:).', 'k-', 'LineWidth', 1.25);
plot(X(:,[1 end]),   Y(:,[1 end]),   'k-', 'LineWidth', 1.25);

% row lines (t = const)
for j = 1:J
    plot(X(j,:), Y(j,:), '-', 'Color', col1, 'LineWidth', 0.8);
end

% column lines (s = const)
for i = 1:I
    plot(X(:,i), Y(:,i), '-', 'Color', col2, 'LineWidth', 0.8);
end

title(sprintf('Vanilla horseshoe/C-grid  Ns=%d, Nt=%d', Ns, Nt));
xlabel('x'); 
ylabel('y');
hold off;

%% ---------------- Optional Jacobian check ----------------
%{
xs = (X(:,3:end)-X(:,1:end-2))/(2*hs);
ys = (Y(:,3:end)-Y(:,1:end-2))/(2*hs);
xtd = (X(3:end,:)-X(1:end-2,:))/(2*ht);
ytd = (Y(3:end,:)-Y(1:end-2,:))/(2*ht);

Jc = xs(2:end-1,:).*ytd(:,2:end-1) - xtd(:,2:end-1).*ys(2:end-1,:);
fprintf('Min det(J) over cells: %.3e\n', min(Jc(:)));
%}
end
