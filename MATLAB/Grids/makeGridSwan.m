function [X, Y] = makeGridSwan(nptsx, nptsy)
% makeGridSwan
% Builds a structured grid inside a "swan-like" domain
% using the boundary curves:
%   f1(s) = bottom
%   f2(s) = right
%   f3(s) = top
%   f4(s) = left
%
% Inputs:
%   nptsx - number of cells in xi-direction
%   nptsy - number of cells in eta-direction
%
% Outputs:
%   X, Y  - grid point coordinates

Ni = nptsx;
Nj = nptsy;

xi  = linspace(0,1,Ni+1);
eta = linspace(0,1,Nj+1);

[X, Y] = meshgrid(xi, eta);

%% Boundary curves
f1 = @(s) [s, 0];                         % bottom
f2 = @(s) [1 + 2*s - 2*s.^2, s];         % right
f3 = @(s) [s, 1 - 3*s + 3*s.^2];         % top
f4 = @(s) [0, s];                        % left

%% Corner points
Pbl = f1(0);   % bottom-left
Pbr = f1(1);   % bottom-right
Ptl = f3(0);   % top-left
Ptr = f3(1);   % top-right

%% Transfinite interpolation
for j = 1:Nj+1
    t = eta(j);

    for i = 1:Ni+1
        s = xi(i);

        % Boundary values
        B = f1(s);   % bottom
        R = f2(t);   % right
        T = f3(s);   % top
        L = f4(t);   % left

        % TFI blending
        P = (1-t)*B + t*T + (1-s)*L + s*R ...
          - (1-s)*(1-t)*Pbl ...
          - s*(1-t)*Pbr ...
          - (1-s)*t*Ptl ...
          - s*t*Ptr;

        X(j,i) = P(1);
        Y(j,i) = P(2);
    end
end

%% Plot the grid
figure('Color','w');
hold on

% Horizontal grid lines
for j = 1:Nj+1
    plot(X(j,:), Y(j,:), 'k-', 'LineWidth', 0.8);
end

% Vertical grid lines
for i = 1:Ni+1
    plot(X(:,i), Y(:,i), 'k-', 'LineWidth', 0.8);
end

axis equal
axis off
title('Swan Shaped Structured Grid', 'FontSize', 14);