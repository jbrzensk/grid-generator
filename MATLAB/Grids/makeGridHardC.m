function [ X, Y ] = makeGridHardC( nptsx, nptsy )
%% Sharp-cornered C boundary in MATLAB

% --- Parameters (example values) ---
a1 = 0.5; b1 = 0.05;  % bottom corner sizes
a2 = 1.0; b2 = 0.50;  % top corner sizes
m = nptsx;            % points along bottom/top boundary
n = nptsy;            % points along left/right boundary

% --- Normalized parameters (0 to 1) ---
tb = linspace(0,1,m+1);   % bottom boundary
tt = linspace(0,1,m+1);   % top boundary
tl = linspace(0,1,n+1);   % left/right boundary

% --- Compute reference fractions ---
xlb = 2*(a1+b1);
xlt = 2*(a2+b2);
rb1 = a1 / xlb;
rb2 = (a1+2*b1)/xlb;
rt1 = a2 / xlt;
rt2 = (a2+2*b2)/xlt;

% --- Initialize arrays ---
xb_bdry = zeros(size(tb));
yb_bdry = zeros(size(tb));
xt_bdry = zeros(size(tt));
yt_bdry = zeros(size(tt));
xl_bdry = zeros(size(tl));
yl_bdry = zeros(size(tl));
xr_bdry = zeros(size(tl));
yr_bdry = zeros(size(tl));

% --- Bottom and top boundaries (piecewise linear) ---
for i = 1:length(tb)
    % Bottom
    if tb(i) <= rb1
        xb_bdry(i) = xlb*tb(i);
        yb_bdry(i) = b1;
    elseif tb(i) > rb1 && tb(i) < rb2
        xb_bdry(i) = a1;
        yb_bdry(i) = xlb*(0.5-tb(i));
    else
        xb_bdry(i) = xlb*(1 - tb(i));
        yb_bdry(i) = -b1;
    end
    
    % Top
    if tt(i) <= rt1
        xt_bdry(i) = xlt*tt(i);
        yt_bdry(i) = b2;
    elseif tt(i) > rt1 && tt(i) < rt2
        xt_bdry(i) = a2;
        yt_bdry(i) = xlt*(0.5 - tt(i));
    else
        xt_bdry(i) = xlt*(1 - tt(i));
        yt_bdry(i) = -b2;
    end
end

% --- Left and right boundaries ---
for j = 1:length(tl)
    xl_bdry(j) = 0;
    yl_bdry(j) = (1 - tl(j))*b1 + tl(j)*b2;
    xr_bdry(j) = 0;
    yr_bdry(j) = -yl_bdry(j);
end



%% Plot the boundaries for visualization
figure('Color','w'); hold on
plot(xb_bdry, yb_bdry, 'b-', 'LineWidth',1.5, 'DisplayName','Bottom') % bottom
plot(xt_bdry, yt_bdry, 'r-', 'LineWidth',1.5, 'DisplayName','Top') % top
plot(xl_bdry, yl_bdry, 'g-', 'LineWidth',1.5, 'DisplayName','Left') % left
plot(xr_bdry, yr_bdry, 'm-', 'LineWidth',1.5, 'DisplayName','Right') % right
axis equal
axis tight
grid on
legend
title('Sharp-cornered C boundary')

%% Interpolate to find the middles
X = zeros(n+1, m+1);
Y = zeros(n+1, m+1);

for j = 1:n+1
    t = tl(j);  % normalized vertical coordinate
    
    for i = 1:m+1
        s = tb(i); % normalized horizontal coordinate
        
        % Bilinear-like transfinite interpolation
        X(j,i) = (1-s)*xl_bdry(j) + s*xr_bdry(j) + ...
                 (1-t)*xb_bdry(i) + t*xt_bdry(i) - ...
                 ((1-s)*(1-t)*xl_bdry(1) + s*(1-t)*xr_bdry(1) + ...
                  (1-s)*t*xl_bdry(end) + s*t*xr_bdry(end));
              
        Y(j,i) = (1-s)*yl_bdry(j) + s*yr_bdry(j) + ...
                 (1-t)*yb_bdry(i) + t*yt_bdry(i) - ...
                 ((1-s)*(1-t)*yl_bdry(1) + s*(1-t)*yr_bdry(1) + ...
                  (1-s)*t*yl_bdry(end) + s*t*yr_bdry(end));
    end
end

%% --- Plot the grid ---
figure('Color','w'); hold on
for j = 1:n+1
    plot(X(j,:), Y(j,:), 'k-', 'LineWidth',0.7); % horizontal grid lines
end
for i = 1:m+1
    plot(X(:,i), Y(:,i), 'k-', 'LineWidth',0.7); % vertical grid lines
end
axis equal; grid on
title('Structured C-shaped Grid (Interior + Boundaries)', 'FontSize',14)