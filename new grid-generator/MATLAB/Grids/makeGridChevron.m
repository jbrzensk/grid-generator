function [ X, Y ] = makeGridChevron( nptsx, nptsy )
% Makes the C grid

Ni = nptsx;  
Nj = nptsy; 

xi  = linspace(0,1,Ni+1);      % xi logical coord
eta = linspace(0,1,Nj+1);      % eta logical coord

[X, Y] = meshgrid(xi, eta);

X(2:end-1, 2:end-1) = 0;
Y(2:end-1, 2:end-1) = 0;

%% Bottom Chevron
xpt = X(end,:);
ynew = abs(xpt - 0.5) + 0.5;
Y(end,:) = ynew;

%% Top Chevron
xpt = X(1,:);
slope = 1; % change this for more/less feature
ynew = slope * (abs(xpt - 0.5) - 0.5);
Y(1,:) = ynew;

% Quasi-interpolation is close enough for show
for j = 1:Nj+1
    t = eta(j);  % normalized vertical position
    % interpolate along the horizontal direction
    for i = 1:Ni+1
        s = xi(i);        % normalized horizontal position
        % linear interpolation between left and right boundaries at this eta
        X(j,i) = (1 - s) * X((j),1) + s * X((j),end);
        % interpolate Y between bottom and top boundaries at this xi
        Y(j,i) = (1 - t) * Y(1,(i)) + t * Y(end,(i));
    end
end
%% --- Plot the grid ---
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
title('Chevron Shaped Structured Grid ( use even nx )', 'FontSize', 14);