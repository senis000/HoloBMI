%pool_stim_files(files_cell)

addpath('D:\Dropbox\Code\export_fig\export_fig_3.3.20'); 


screensize = get(0, 'ScreenSize')
%1) Load the psth data, decide on stim lag
idx = 0;
file_cell = {}

idx = idx+1; 
file_cell{idx} = 'D:\Dropbox\Data\grid_stim\200310_finegrid\fine_gridstim-000\psth_data.mat';
idx = idx+1; 
file_cell{idx} = 'D:\Dropbox\Data\grid_stim\200310_finegrid\fine_gridstim-001\psth_data.mat';
idx = idx+1; 
file_cell{idx} = 'D:\Dropbox\Data\grid_stim\200310_finegrid\fine_gridstim-002\psth_data.mat';
num_files = length(file_cell)
for i = 1:num_files
    exist(file_cell{i})
end
save_dir = 'D:\Dropbox\Data\grid_stim\200310_finegrid';
pool_data_path = fullfile('D:\Dropbox\Data\grid_stim\200310_finegrid', 'psth_pool_data.mat'); 
grid_data_path = fullfile('D:\DATA\grid_stim\200310_finegrid\NVI20\D0', 'fine_grid.mat'); 

%%
%Plot an individual file's psth with an applied shift.
%Re plot the analysis data with the applied shift.
%Pool data across files. 
%%
%Data2save: 
%target_psth_mat
%target_psth_x
%grid_psth_mat
%grid_psth_x
data2pool = repmat(struct(...
    'lag',[], ...
    'target_psth_mat',[], ...
    'target_psth_x',[], ...
    'over_grid',[], ...
    'grid_psth_x',[]), [num_files 1])

%%
clear dff_psth
clear psth_plot
clear response_over_grid

idx_load = 1; 
load(file_cell{idx_load}); 

%idx=1, lag=4
%idx=2, lag=3
%idx=3, lag=4
%%
test = load(file_cell{idx_load})
%%
%Lagged PSTH: 
%--------------------------------------------------------------------------
stim_lag = 4
%--------------------------------------------------------------------------
psth_plot = dff_psth; 

zoom_win_t = [-4.5 4.5]
zoom_win = zoom_win_t*30

zoom_win_o = zoom_win-stim_lag
zoom_win_o_t = zoom_win_o/30
sel = find((psth_plot.x >= zoom_win_o_t(1)) & (psth_plot.x <= zoom_win_o_t(2)));

%DATA
target_psth_mat     = psth_plot.mat(sel,:,:);
target_psth_x       = (zoom_win(1):zoom_win(2))/30; 
%ASSIGN:
data2pool(idx_load).lag                 = stim_lag; 
data2pool(idx_load).target_psth_mat     = target_psth_mat; 
data2pool(idx_load).target_psth_x       = target_psth_x; 

%PLOT
target_psth_mean    = psth_plot.mean(sel);
h = figure;
hold on;
for i =1:size(target_psth_mat,3)
    plot(target_psth_x, target_psth_mat(:,:,i)); 
end
set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('dff'); 
title('dff locked to stim at target (stim at t=0)')
vline(0)
plot(target_psth_x, target_psth_mean, 'k', 'LineWidth', 1.3)

%%
%--------------------------------------------------------------------------
%Grid PSTH

psth_win    = grid_psth_win; 
x_plot      = (psth_win(1):psth_win(2))/30; 
sel         = find((x_plot >= zoom_win_o_t(1)) & (x_plot <= zoom_win_o_t(2)));

%DATA
grid_psth_x     = (zoom_win(1):zoom_win(2))/30; 
over_grid       = response_over_grid; 
%Select idxs from the grid data: 
for grid_idx = 1:num_grid_pts
    over_grid(grid_idx).psth_mat = ...
        over_grid(grid_idx).psth_mat(sel,:,:);
    over_grid(grid_idx).psth_mean = ...
        over_grid(grid_idx).psth_mean(sel,:);
end
%ASSIGN:
data2pool(idx_load).over_grid           = over_grid;
data2pool(idx_load).grid_psth_x         = grid_psth_x;

incr_width = 1/(num_x_incr+1); 
incr_height = 1/(num_y_incr+1); 
%Range can be changed: 
ymin = -0.1;
ymax = 1.5; 
plot_individual_trials = 1; 
h = figure('Position', [screensize(3)/2 1 screensize(3)/4 screensize(3)/4]);%divide by 2.5 to have gray space
% h = figure; 
for grid_idx = 1:num_grid_pts %num_grid_pts%1:num_grid_pts
    ax = axes('Position', [x_grid_norm(grid_idx) y_grid_norm(grid_idx) incr_width, incr_height]); 
    hold on; 
    plot([0 0], [ymin ymax], 'r'); 
    
    if plot_individual_trials
        mat_plot = over_grid(grid_idx).psth_mat; 
        num_trials = size(mat_plot,3); 
        for j = 1:num_trials
            plot(grid_psth_x, mat_plot(:,1,j), 'Color', 0.8*ones(3,1)); 
        end
        plot(grid_psth_x, over_grid(grid_idx).psth_mean, 'k'); 
    end    
    ylim([ymin ymax]); 
    xlim(zoom_win_t)
    set(gca, 'xticklabel',[])
    set(gca, 'yticklabel',[])
%     axis square
end

%%
%Plot all PSTH
for idx = 1:num_files 
    x_plot      = data2pool(idx).target_psth_x; 
    mat         = data2pool(idx).target_psth_mat; 
    h = figure;
    hold on;
    for i =1:size(mat,3)
        plot(x_plot, mat(:,1,i)); 
    end
    set(gca,'TickDir','out');
    xlabel('time (s)'); 
    ylabel('dff'); 
    title('dff locked to stim at target (stim at t=0)')
    vline(0)
    plot(x_plot, mean(mat,3), 'k', 'LineWidth', 1.3)
end

%%
%Plot grid: 
incr_width = 1/(num_x_incr+1); 
incr_height = 1/(num_y_incr+1); 
%Range can be changed: 
ymin = -0.1;
ymax = 1.5; 
plot_individual_trials = 1;

for idx = 1:num_files
    grid_psth_x     = data2pool(idx).grid_psth_x; 
    over_grid       = data2pool(idx).over_grid; 

    h = figure('Position', [screensize(3)/2 1 screensize(3)/4 screensize(3)/4]);%divide by 2.5 to have gray space
    % h = figure; 
    for grid_idx = 1:num_grid_pts %num_grid_pts%1:num_grid_pts
        ax = axes('Position', [x_grid_norm(grid_idx) y_grid_norm(grid_idx) incr_width, incr_height]); 
        hold on; 
        plot([0 0], [ymin ymax], 'r'); 
        if plot_individual_trials
            mat_plot = over_grid(grid_idx).psth_mat; 
            num_trials = size(mat_plot,3); 
            for j = 1:num_trials
                plot(grid_psth_x, mat_plot(:,1,j), 'Color', 0.8*ones(3,1)); 
            end
            plot(grid_psth_x, over_grid(grid_idx).psth_mean, 'k'); 
        end    
        ylim([ymin ymax]); 
        xlim(zoom_win_t)
        set(gca, 'xticklabel',[])
        set(gca, 'yticklabel',[])
    %     axis square
    end
end

%%
%Pool data over files
save_bool = 0
if save_bool
    data2pool_path = fullfile(save_dir, 'psth_data2pool.mat'); 
    save(data2pool_path, 'data2pool'); 
end

%%
%Pool the data over files: 
pool_data = data2pool(1); 
for idx = 2:num_files
    pool_data.target_psth_mat = ...
        cat(3, pool_data.target_psth_mat, data2pool(idx).target_psth_mat);
    for grid_idx = 1:num_grid_pts
        pool_data.over_grid(grid_idx).psth_mat = ...
            cat(3, ...
            pool_data.over_grid(grid_idx).psth_mat, ...
            data2pool(idx).over_grid(grid_idx).psth_mat); 
    end   
end
for grid_idx = 1:num_grid_pts
    pool_data.over_grid(grid_idx).psth_mean = ...
        mean(pool_data.over_grid(grid_idx).psth_mat,3);
end

save_bool = 0
if save_bool
    pool_data_path = fullfile(save_dir, 'psth_pool_data.mat'); 
    save(pool_data_path, 'pool_data'); 
end

%%
%Load pool data and grid data
load(grid_data_path)
load(pool_data_path)

x_step_vec = x_step:x_step:x_extent(2);
y_step_vec = y_step:y_step:y_extent(2);

num_x_incr = length(x_step_vec)*2; 
num_y_incr = length(y_step_vec)*2; 

x_grid_flat = x_mesh_flat; 
y_grid_flat = (dim - y_mesh_flat + 1); 

x_min   = min(x_grid_flat);
x_range = max(x_grid_flat)-min(x_grid_flat); 
y_range = max(y_grid_flat)-min(y_grid_flat); 
y_min   = min(y_grid_flat); 

x_grid_norm = (x_grid_flat-x_min)/((num_x_incr+1)/num_x_incr*x_range); 
y_grid_norm = (y_grid_flat-y_min)/((num_y_incr+1)/num_y_incr*y_range); 

incr_width = 1/(num_x_incr+1); 
incr_height = 1/(num_y_incr+1); 

%%
%PSTH: 
save_bool   = 0;

x_plot      = pool_data.target_psth_x; 
mat         = pool_data.target_psth_mat; 
h = figure;
hold on;
for i =1:size(mat,3)
    plot(x_plot, mat(:,1,i)); 
end
axis square
set(gca,'TickDir','out');
% xlabel('time (s)'); 
% ylabel('dff'); 
% title('dff locked to stim at target (stim at t=0)')
vline(0)
plot(x_plot, mean(mat,3), 'k', 'LineWidth', 1.3)

if save_bool
    save_path = fullfile(save_dir, 'target_psth.png');
    export_fig(save_path); 
    save_path = fullfile(save_dir, 'target_psth.eps');
    export_fig(save_path);
end
%%
%GRID PSTH: 
save_bool = 0;
grid_psth_x     = pool_data.grid_psth_x; 
over_grid       = pool_data.over_grid; 

ymin = -0.1;
ymax = 1.5; 
plot_individual_trials = 1; 
h = figure('Position', [screensize(3)/2 1 screensize(3)/3 screensize(3)/3]);%divide by 2.5 to have gray space
% h = figure; 
for grid_idx = 1:num_grid_pts %num_grid_pts%1:num_grid_pts
    ax = axes('Position', [x_grid_norm(grid_idx) y_grid_norm(grid_idx) incr_width, incr_height]); 
    hold on; 
    plot([0 0], [ymin ymax], 'r'); 
    if plot_individual_trials
        mat_plot = over_grid(grid_idx).psth_mat; 
        num_trials = size(mat_plot,3); 
        for j = 1:num_trials
            plot(grid_psth_x, mat_plot(:,1,j), 'Color', 0.8*ones(3,1)); 
        end
        plot(grid_psth_x, over_grid(grid_idx).psth_mean, 'k'); 
    end    
    ylim([ymin ymax]); 
    xlim([min(grid_psth_x) max(grid_psth_x)])
    set(gca, 'xticklabel',[])
    set(gca, 'yticklabel',[])
%     axis square
end
if save_bool
    save_path = fullfile(save_dir, 'grid_psth.png');
    export_fig(save_path); 
    save_path = fullfile(save_dir, 'grid_psth.eps');
    export_fig(save_path); 
end

%%
%Re analyze and plot results:

%%
%Recalculate change over grid: 
change_over_grid = repmat(struct(...
    'pre_mat', [], ...
    'post_mat', [], ...
    'perc_change_mat', [], ...
    'perc_change', []), [num_grid_pts, 1]); 

pre_t = -0.5
post_t = 0.5

x_plot      = grid_psth_x; 
pre_idxs    = find((x_plot >= pre_t & x_plot < 0)); 
post_idxs   = find((x_plot >= 0 & x_plot <= post_t));

%
for grid_idx = 1:num_grid_pts
    mat_i = over_grid(grid_idx).psth_mat; 
    pre_mat_i = squeeze(mean(mat_i(pre_idxs, 1, :), 1)); 
    post_mat_i = squeeze(mean(mat_i(post_idxs, 1, :), 1)); 
    perc_mat_i = (post_mat_i-pre_mat_i)./pre_mat_i; 
    perc_i = 100*mean(perc_mat_i); 
    perc_mean_i = 100*(mean(post_mat_i)-mean(pre_mat_i))/mean(pre_mat_i);
    delta_mean_i = 100*(mean(post_mat_i)-mean(pre_mat_i));
    
    change_over_grid(grid_idx).pre_mat = pre_mat_i; 
    change_over_grid(grid_idx).post_mat = post_mat_i; 
    change_over_grid(grid_idx).perc_change_mat = perc_mat_i; 
    change_over_grid(grid_idx).perc_change = perc_i;
    change_over_grid(grid_idx).perc_change_mean = perc_mean_i;
    change_over_grid(grid_idx).delta_mean = delta_mean_i;
end

%
%Response over Distance: 
%Make a vector of stim distance to target point, and the percent change of
%response, trial averaged
d_px    = []; 
d_um    = [];
r_delta = []; 
r_perc  = [];
%Target: x = col, y = row
for grid_idx = 1:num_grid_pts
    d_px_i = sqrt((over_grid(grid_idx).col - target_roi.x)^2 + (over_grid(grid_idx).row - target_roi.y)^2);
    d_um_i = sqrt(micronsPerPixel.x^2*(over_grid(grid_idx).col - target_roi.x)^2 + micronsPerPixel.y^2*(over_grid(grid_idx).row - target_roi.y)^2);
    d_px = [d_px d_px_i]; 
    d_um = [d_um d_um_i]; 
    
    r_delta_i   = change_over_grid(grid_idx).delta_mean;
    r_delta     = [r_delta r_delta_i]; 
    
    r_perc_i    = change_over_grid(grid_idx).perc_change_mean;
    r_perc      = [r_perc r_perc_i]; 
end

%%
%FOV with square drawn on grid region
%Plot a square around the grid
%xmin: 87
%xmax: 143
%ymin: 231
%ymax: 287
border = 7;
sq_xmin = min(x_mesh_flat)-border; 
sq_xmax = max(x_mesh_flat)+border; 
sq_ymin = min(y_mesh_flat)-border; 
sq_ymax = max(y_mesh_flat)+border; 

bg_sq = target_roi_data.im_bg; 
%Rows: 
bg_sq(sq_ymin:sq_ymax, sq_xmin, :) = 1; 
bg_sq(sq_ymin:sq_ymax, sq_xmax, :) = 1; 
%Cols: 
bg_sq(sq_ymin, sq_xmin:sq_xmax, :) = 1; 
bg_sq(sq_ymax, sq_xmin:sq_xmax, :) = 1; 

h = figure; 
imagesc(bg_sq);
axis square
hold on; 
% scatter(target_roi.x, target_roi.y+3, 700, 'w')
set(h, 'Position', screensize); 
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
save_path = fullfile(save_dir, 'FOV_grid_sq.png');
export_fig(save_path); 
save_path = fullfile(save_dir, 'FOV_grid_sq.eps');
export_fig(save_path); 
% for grid_idx = 1:num_grid_pts
%     val = change_over_grid(grid_idx).delta_mean;
%     scatter(x_mesh_flat(grid_idx), y_mesh_flat(grid_idx), 20, val, 'filled') %, 'x', 'LineWidth', 1.5); 
% end

%%
%Plot zoom in: 
save_bool = 1; 

bg_zoom = target_roi_data.im_bg(sq_ymin:sq_ymax, sq_xmin:sq_xmax,:); 
[bg_zoom_sc, zoom_min, zoom_max] = scale_im(bg_zoom, 0, 98.);
h = figure;
imagesc(bg_zoom_sc); 
axis square
hold on; 
for grid_idx = 1:num_grid_pts
    val = change_over_grid(grid_idx).delta_mean;
    scatter(x_mesh_flat(grid_idx)-sq_xmin, y_mesh_flat(grid_idx)-sq_ymin, 100, val, 'filled') %, 'x', 'LineWidth', 1.5); 
end
colormap('jet'); 
colorbar; 

set(h, 'Position', screensize); 
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

if save_bool
    save_path = fullfile(save_dir, 'grid_response_change.png');
    export_fig(save_path); 
    save_path = fullfile(save_dir, 'grid_response_change.eps');
    export_fig(save_path); 
end

%%
%Plot zoom in: 
save_bool = 1; 

bg_zoom = target_roi_data.im_bg(sq_ymin:sq_ymax, sq_xmin:sq_xmax,:); 
[bg_zoom_sc, zoom_min, zoom_max] = scale_im(bg_zoom, 0, 98.);
h = figure;
hold on; 
for grid_idx = 1:num_grid_pts
    val = change_over_grid(grid_idx).delta_mean;
    scatter(x_mesh_flat(grid_idx)-sq_xmin, y_mesh_flat(grid_idx)-sq_ymin, 100, val, 'filled') %, 'x', 'LineWidth', 1.5); 
end
axis square
colormap('jet'); 
colorbar; 

% set(h, 'Position', screensize); 
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

if save_bool
    save_path = fullfile(save_dir, 'grid_response_change_no_bg.png');
    export_fig(save_path); 
    save_path = fullfile(save_dir, 'grid_response_change_no_bg.eps');
    export_fig(save_path); 
end

%%
h = figure; 
imagesc(target_roi_data.im_bg);
hold on;
scatter(target_roi.x, target_roi.y+3, 600, 'w')
axis square 
set(h, 'Position', screensize); 
%%
%change in dff: 
%--------------------------------------------------------------------------
scatter_size = 400; 
h = figure; 
imagesc(target_roi_data.im_bg);
hold on; 
for grid_idx = 1:num_grid_pts
    val = change_over_grid(grid_idx).delta_mean;
    scatter(x_mesh_flat(grid_idx), y_mesh_flat(grid_idx), 20, val, 'filled') %, 'x', 'LineWidth', 1.5); 
end
colormap('jet'); 
% caxis([0 2000]); 
colorbar
axis square
set(h, 'Position', screensize); 

% title('Percent change in DFF of target ROI for stims marked on image'); 
% if save_bool
%     filename = 'stim_response_on_bg_perc_change'
%     save_ext(save_dir, filename, ext_cell)
% end

%%
h= figure; 
imagesc(target_roi_data.im_roi)
axis square
%%
%--------------------------------------------------------------------------
%Plot change in DFF over the grid points: 

%%
d = d_um;
r = r_delta;

d_vec = unique(d); 
r_vec = zeros(size(d_vec)); 
for i =1:length(d_vec)
    sel = find(d==d_vec(i)); 
    r_vec(i) = mean(r(sel)); 
end

save_bool = 1
h = figure;
hold on
scatter(d,r, 'fill'); 
plot(d_vec, r_vec, 'k')
xlim([-5 35]); 
xlabel('distance between grid stim point and target roi (um)'); 
ylabel('change in DFF, 0.5s before and after stim'); 
axis square
set(gca,'TickDir','out');

if save_bool
    plot_name = 'response_vs_distance'
    tic
    disp('saving...')
    save_path = fullfile(save_dir, [plot_name '.eps']);
%     saveas(h, save_path); 
    export_fig(save_path); 
    
    save_path = fullfile(save_dir, [plot_name '.png']);
%     saveas(h, save_path); 
    export_fig(save_path);
    disp('DONE saving'); 
    toc
else
    disp('not saving...'); 
end

%%
h = figure;
scatter(d, r, 'fill'); 


%%
h = figure;
hold on; 
scatter(d, r, 'fill'); 
plot(d_vec, r_vec, 'k')
xlabel('distance between grid stim point and target roi (um)'); 
ylabel('change in DFF, 0.5s before and after stim'); 
axis square
set(gca,'TickDir','out');


