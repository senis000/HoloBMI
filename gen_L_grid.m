function [x_mesh_flat, y_mesh_flat, ctr_idx, num_grid_pts] = ...
    gen_L_grid(target, dim, x_extent, y_extent, x_step, y_step)

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

remove_x = find(...
    (target.x+x_step_vec) < x_extent_centered(1) | ...
    (target.x+x_step_vec) > x_extent_centered(2))
x_step_vec(remove_x) = [];

remove_y = find(...
    (target.y+y_step_vec) < y_extent_centered(1) | ...
    (target.y+y_step_vec) > y_extent_centered(2))
y_step_vec(remove_y) = [];

x_mesh_flat = target.x + [0 -x_step_vec zeros(size(y_step_vec))]; 
y_mesh_flat = target.y + [0 zeros(size(x_step_vec)) -y_step_vec];
 
ctr_idx = find(x_mesh_flat == target.x & y_mesh_flat == target.y); 
num_grid_pts = length(x_mesh_flat);
