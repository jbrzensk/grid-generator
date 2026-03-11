function [ X, Y ] = makeGridC( nptsx, nptsy )
% Makes the C grid

Ni = nptsx;  
Nj = nptsy; 

xi  = linspace(0,1,Ni+1);      % xi logical coord
eta = linspace(0,1,Nj+1);      % eta logical coord

[Xi, Eta] = meshgrid(xi, eta);

% --- Parameters (from Fortran code) ---
pi_val = pi;
rho = 2.0;
b0 = 1.0;
b1 = 2.0;

% --- Compute physical coordinates ---
R = b0 + (b1 - b0) * Eta;            % radial distance
Theta = pi_val * (1 - 2*Xi)/2;       % angular coordinate

X = rho * R .* cos(Theta);           % x-coordinate
Y = R .* sin(Theta);                  % y-coordinate

% --- Plot the grid ---
figure('Color','w'); hold on

% Radial lines (eta-direction)
for j = 1:Nj+1
    plot(X(j,:), Y(j,:), 'k-', 'LineWidth', 0.8);
end

% Circular lines (xi-direction)
for i = 1:Ni+1
    plot(X(:,i), Y(:,i), 'k-', 'LineWidth', 0.8);
end

axis equal
axis off
title('Horseshoe (C-shaped) Structured Grid', 'FontSize', 14);