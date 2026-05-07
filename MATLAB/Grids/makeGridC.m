function [ X, Y ] = makeGridC( nptsx, nptsy )
% Makes the smooth / soft C grid and plots its boundary curves

Ni = nptsx;
Nj = nptsy;

xi  = linspace(0,1,Ni+1);
eta = linspace(0,1,Nj+1);

[Xi, Eta] = meshgrid(xi, eta);

%% --- Parameters ---
rho = 2.0;
b0 = 1.0;
b1 = 2.0;

%% --- Compute physical coordinates ---
R = b0 + (b1 - b0) * Eta;
Theta = pi * (1 - 2*Xi)/2;

X = rho * R .* cos(Theta);
Y = R .* sin(Theta);

%% --- Extract boundaries ---
xb_bdry = X(1,:);
yb_bdry = Y(1,:);

xt_bdry = X(end,:);
yt_bdry = Y(end,:);

xl_bdry = X(:,1);
yl_bdry = Y(:,1);

xr_bdry = X(:,end);
yr_bdry = Y(:,end);

%% --- Plot boundary only ---
figure('Color','w'); hold on

plot(xb_bdry, yb_bdry, 'b-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Bottom / inner arc');

plot(xt_bdry, yt_bdry, 'r-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Top / outer arc');

plot(xl_bdry, yl_bdry, 'g-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Left cut');

plot(xr_bdry, yr_bdry, 'm-', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Right cut');

axis equal
axis tight
grid on
legend
title('Smooth C boundary')

%% --- Plot the full grid ---
figure('Color','w'); hold on

for j = 1:Nj+1
    plot(X(j,:), Y(j,:), 'k-', 'LineWidth', 0.8);
end

for i = 1:Ni+1
    plot(X(:,i), Y(:,i), 'k-', 'LineWidth', 0.8);
end

axis equal
grid on
title('Smooth C-shaped Structured Grid', 'FontSize', 14);

end