function [ X, Y ] = makeGridFish( nptsx, nptsy )
% Makes the C grid

Ni = nptsx;  
Nj = nptsy; 

xi  = linspace(0,1,Ni+1);      % xi logical coord
eta = linspace(0,1,Nj+1);      % eta logical coord

[X, Y] = meshgrid(xi, eta);

X(2:end-1, 2:end-1) = 0;
Y(2:end-1, 2:end-1) = 0;

%% Right Bulge in
xpt = X(:,end);
ypt = Y(:,end);

amp = -0.1;
xnew = amp .* cos(2*pi .* ypt) - amp;
xnew = xpt - xnew;
X(:,end) = xnew;

%% Bottom Dome
xpt = X(1,:);
amp = -0.2;  % Change this for more/less feature
ynew = amp .* cos(2*pi .* xpt)  - amp;
Y(1,:) = ynew;

%% Top Valley
xpt = X(end,:);
amp = 0.2; % Change this for more/less feature
ynew = ( amp .* cos(2*pi .* xpt) ) + (1 - amp);
Y(end,:) = ynew;

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
title('Fish Shaped Structured Grid', 'FontSize', 14);