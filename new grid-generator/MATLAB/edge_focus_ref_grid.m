%% Edge focus reference grid generation
%
% BY: Jared Brzenski
%
%
%
cluster_factor = 5; % Higher values result in tighter clustering 
grid_range_min = 0;
grid_range_max = 1;

grid_x = makeAsinhVec( nx, cluster_factor, grid_range_min, grid_range_max);
grid_y = makeAsinhVec( ny, cluster_factor, grid_range_min, grid_range_max);

% --- Create the 2D grid using meshgrid ---
[Xl, Yl] = meshgrid(grid_x, grid_y);

% --- Visualize the grid ---
figure;
scatter(Xl(:), Yl(:), 'filled');
title('2D Reference Grid');
xlabel('x');
ylabel('y');
grid on;
axis equal;
xlim([0 1]);
ylim([0 1]);


function [grid_x] = makeAsinhVec( num_points, cluster_factor, grid_range_min, grid_range_max)
    % --- Create the 1D clustered vector using asinh ---
    linear_points = linspace(-1, 1, num_points);
    asinh_points = asinh(linear_points * cluster_factor);
    
    % --- Manually scale the points ---
    input_min = min(asinh_points);
    input_max = max(asinh_points);
    output_min = grid_range_min;
    output_max = grid_range_max;
    
    if input_max == input_min
        % Handle the constant input case (i.e., all points are the same)
        % This can happen if cluster_factor is very small.
        grid_x = ones(1, num_points) * (output_min + output_max) / 2;
    else
        % Perform the scaling manually
        grid_x = output_min + (asinh_points - input_min) * (output_max - output_min) / (input_max - input_min);
    end
end