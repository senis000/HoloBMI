%pool_stim_files(files_cell)
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

idx_load = 3; 
load(file_cell{idx_load}); 

%idx=1, lag=4
%idx=2, lag=3
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

save_bool = 1
if save_bool
    pool_data_path = fullfile(save_dir, 'psth_pool_data.mat'); 
    save(pool_data_path, 'pool_data'); 
end
%%
%PSTH: 
x_plot      = pool_data.target_psth_x; 
mat         = pool_data.target_psth_mat; 
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


%GRID PSTH: 
grid_psth_x     = pool_data.grid_psth_x; 
over_grid       = pool_data.over_grid; 

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
%Re analyze and plot results:

% %%
% %Recalculate change over grid: 
% change_over_grid = repmat(struct(...
%     'pre_mat', [], ...
%     'post_mat', [], ...
%     'perc_change_mat', [], ...
%     'perc_change', []), [num_grid_pts, 1]); 
% 
% pre_idxs    = find((x_plot >= -2 & x_plot < 0)); 
% post_idxs   = find((x_plot >= 0 & x_plot <= 2));
% 
% %
% for grid_idx = 1:num_grid_pts
%     mat_i = response_over_grid(grid_idx).psth_mat; 
%     pre_mat_i = squeeze(mean(mat_i(pre_idxs, 1, :), 1)); 
%     post_mat_i = squeeze(mean(mat_i(post_idxs, 1, :), 1)); 
%     perc_mat_i = (post_mat_i-pre_mat_i)./pre_mat_i; 
%     perc_i = 100*mean(perc_mat_i); 
%     perc_mean_i = 100*(mean(post_mat_i)-mean(pre_mat_i))/mean(pre_mat_i);
%     delta_mean_i = 100*(mean(post_mat_i)-mean(pre_mat_i));
%     
%     change_over_grid(grid_idx).pre_mat = pre_mat_i; 
%     change_over_grid(grid_idx).post_mat = post_mat_i; 
%     change_over_grid(grid_idx).perc_change_mat = perc_mat_i; 
%     change_over_grid(grid_idx).perc_change = perc_i;
%     change_over_grid(grid_idx).perc_change_mean = perc_mean_i;
%     change_over_grid(grid_idx).delta_mean = delta_mean_i;
% end
% 
% %%
% %Response over Distance: 
% %Make a vector of stim distance to target point, and the percent change of
% %response, trial averaged
% d_px    = []; 
% d_um    = [];
% r_delta = []; 
% r_perc  = [];
% %Target: x = col, y = row
% for grid_idx = 1:num_grid_pts
%     d_px_i = sqrt((response_over_grid(grid_idx).col - target_roi.x)^2 + (response_over_grid(grid_idx).row - target_roi.y)^2);
%     d_um_i = sqrt(micronsPerPixel.x^2*(response_over_grid(grid_idx).col - target_roi.x)^2 + micronsPerPixel.y^2*(response_over_grid(grid_idx).row - target_roi.y)^2);
%     d_px = [d_px d_px_i]; 
%     d_um = [d_um d_um_i]; 
%     
%     r_delta_i   = change_over_grid(grid_idx).delta_mean;
%     r_delta     = [r_delta r_delta_i]; 
%     
%     r_perc_i    = change_over_grid(grid_idx).perc_change_mean;
%     r_perc      = [r_perc r_perc_i]; 
% end

