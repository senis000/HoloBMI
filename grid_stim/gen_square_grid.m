function [x_mesh_flat, y_mesh_flat, ctr_idx, num_grid_pts, x_mesh, y_mesh] = ...
    gen_square_grid(target, dim, x_extent, y_extent, x_step, y_step)
%INPUT:
%target - struct with fields x,y.  Grid is centered on the target
%dim - image is dim X dim pixels
% x_extent - 1x2, boundaries of the grid around the target.x 
% y_extent - 1x2, boundaries of the grid around target.y
% x_step -
% y_step - 

x_extent_centered = target.x + x_extent;
x_extent_centered(x_extent_centered < 1) = 1;
x_extent_centered(x_extent_centered > dim) = dim;

y_extent_centered = target.y + y_extent;
y_extent_centered(y_extent_centered < 1) = 1;
y_extent_centered(y_extent_centered > dim) = dim;
%
x_step_vec = x_step:x_step:x_extent(2);
y_step_vec = y_step:y_step:y_extent(2);
%
x_grid = target.x + [sort(-x_step_vec) 0 x_step_vec]; 
remove_x = find(...
    x_grid < x_extent_centered(1) | ...
    x_grid > x_extent_centered(2))
x_grid(remove_x) = []

y_grid = target.y + [sort(-y_step_vec) 0 y_step_vec]; 
remove_y = find(...
    y_grid < y_extent_centered(1) | ...
    y_grid > y_extent_centered(2))
y_grid(remove_y) = []
%Sort is just used for ease of plotting the results
%
[x_mesh, y_mesh] = meshgrid(x_grid, y_grid);
%
x_mesh_flat = x_mesh(:); 
y_mesh_flat = y_mesh(:); 
ctr_idx = find(x_mesh_flat == target.x & y_mesh_flat == target.y); 
num_grid_pts = length(x_mesh_flat);