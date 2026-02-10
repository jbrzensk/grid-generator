function [ X, Y ] = makeGridA( nptsx, nptsy )
% This makes the crazy Gaussian Humped grid 

% Grid size ---
Ni = nptsx;   % number of points in x-direction
Nj = nptsy;   % number of points in y-direction

xi  = linspace(0,1,Ni+1);      % computational coord along x
eta = linspace(0,1,Nj+1);      % computational coord along y

% uncomment below for edge clustering
%xi = (sin(linspace(0,pi/2,Ni))).^2;
%eta = (sin(linspace(0,pi/2,Nj))).^2;

% physical domain base size ---
Lx = 12;
Ly = 8;

% gaussian hump function 
gauss = @(s, pos, wid, amp) amp * exp(-((s - pos).^2) / (2*wid^2));

% Define humps for each side
% Each side: array of [position_normalized, width (in normalized coord), amplitude_physical]
% amplitude: positive = pushes INTO domain (inward normal)

% Bottom (y=0) inward is +y
bottom_humps = [
    0.25, 0.05, 1.0;
    0.65, 0.10, 2.0
];

% Top (y=Ly) inward is -y
top_humps = [
    0.30, 0.08, -1.0;
    0.65, 0.05, -3.6
];

% Left (x=0) inward is +x
left_humps = [
    0.50, 0.10, 3.0
];

% Right (x=Lx) inward is -x
right_humps = [
    0.6, 0.09, -1.1
];

%% No need to edit below here

% --- Build boundary curves as vector functions of normalized parameter
% Bottom boundary (xi -> x_b(xi), y_b(xi))
xb = Lx * xi;                                   % base x along bottom
yb = zeros(size(xi));
for k=1:size(bottom_humps,1)
    yb = yb + gauss(xi, bottom_humps(k,1), bottom_humps(k,2), bottom_humps(k,3));
end

% Top boundary (xi -> x_t(xi), y_t(xi))
xt = Lx * xi;
yt = Ly * ones(size(xi));
for k=1:size(top_humps,1)
    yt = yt + gauss(xi, top_humps(k,1), top_humps(k,2), top_humps(k,3));
end

% Left boundary (eta -> x_l(eta), y_l(eta))
yl = Ly * eta;                % base y along left
xl = zeros(size(eta));
for k=1:size(left_humps,1)
    xl = xl + gauss(eta, left_humps(k,1), left_humps(k,2), left_humps(k,3));
end

% Right boundary (eta -> x_r(eta), y_r(eta))
yr = Ly * eta;
xr = Lx * ones(size(eta));
for k=1:size(right_humps,1)
    xr = xr + gauss(eta, right_humps(k,1), right_humps(k,2), right_humps(k,3));
end

% --- Corner points
P00 = [xb(1), yb(1)];   % bottom-left
P10 = [xb(end), yb(end)]; % bottom-right
P01 = [xt(1), yt(1)];   % top-left
P11 = [xt(end), yt(end)]; % top-right

[X, Y] = deal(zeros(Nj+1, Ni+1));

% Quasi-interpolation is close enough for show
for j = 1:Nj+1
    t = eta(j);  % normalized vertical position
    % interpolate along the horizontal direction
    for i = 1:Ni+1
        s = xi(i);        % normalized horizontal position
        % linear interpolation between left and right boundaries at this eta
        X(j,i) = (1 - s) * xl(j) + s * xr(j);
        % interpolate Y between bottom and top boundaries at this xi
        Y(j,i) = (1 - t) * yb(i) + t * yt(i);
    end
end
% --- Plot results ---
figure('Color','w','Position',[100 100 800 500])
hold on
% draw grid lines
for j = 1:Nj
    plot(X(j,:), Y(j,:), 'k-', 'LineWidth', 0.7);
end
for i = 1:Ni
    plot(X(:,i), Y(:,i), 'b-', 'LineWidth', 0.7);
end

% plot boundary curves thicker
plot(xb, yb, 'r-', 'LineWidth', 2)
plot(xt, yt, 'r-', 'LineWidth', 2)
plot(xl, yl, 'r-', 'LineWidth', 2)
plot(xr, yr, 'r-', 'LineWidth', 2)

axis equal tight
axis off
title('Crazy Grid A', 'FontSize', 12)
drawnow
