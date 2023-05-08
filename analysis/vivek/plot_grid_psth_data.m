function plot_grid_psth_data(psth_data_path, save_dir, save_bool)
load(psth_data_path)
close all
screensize = get(0, 'ScreenSize')
ext_cell = {'.png'} %{'.eps', '.png'}
% %%
% %--------------------------------------------------------------------------
% h = figure('Position', [screensize(3)/2 1 screensize(3)/2 screensize(4)]);
% imagesc(im_bg); 
% colormap('gray'); 
% axis square
% title(['Background Image ' 'min prc: ' num2str(im_bg_min_perc) '; max prc: ' num2str(im_bg_max_perc)]); 
% 
% if save_bool
%     filename = 'im_bg'
%     save_ext(save_dir, filename, ext_cell)
% end
% 
% %%
% %--------------------------------------------------------------------------
% %Grid points marked with 'X'
% h = figure('Position', [screensize(3)/2 1 screensize(3)/2 screensize(4)]);
% imagesc(im_bg); 
% colormap('gray')
% hold on
% scatter(x_mesh_flat, y_mesh_flat, 200, 'x', 'LineWidth', 1.5); 
% scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 200, 'r', 'x', 'LineWidth', 1.5); 
% xlim([1 dim])
% ylim([1 dim])
% axis square
% %Reminder, mapping from x,y to image space: 
% %x = column, y = row
% if save_bool
%     filename = 'grid_x'
%     save_ext(save_dir, filename, ext_cell)
% end
% 
% %%
% %--------------------------------------------------------------------------
% %Number of Stims Histogram
% num_stims = length(stim_sequence)
% num_target_stims = length(find(stim_sequence == ctr_idx))
% %Number of stimulations per sites:
% h= figure;
% hist(stim_sequence, 1:num_grid_pts);
% set(gca,'TickDir','out');
% xlabel('idx of grid point'); 
% ylabel('number of stims to grid point'); 
% title(['num stims per grid point, idx of grid pt at target roi: ' num2str(ctr_idx)])
% disp('center point:')
% ctr_idx
% if save_bool
%     filename = 'stim_hist'
%     save_ext(save_dir, filename, ext_cell)
% end
% %%
% %--------------------------------------------------------------------------
% %Order of stimulation: 
% h = figure; hold on; 
% plot(stim_sequence, 'LineWidth', 0.5)
% scatter(1:num_stims, stim_sequence)
% set(gca,'TickDir','out');
% xlabel('stim idx'); 
% ylabel('grid idx stimmed'); 
% title('order of stimulations')
% if save_bool
%     filename = 'stim_seq'
%     save_ext(save_dir, filename, ext_cell)
% end
% 
%%
%--------------------------------------------------------------------------
%DFF Target ROI PSTH + Individual Traces
psth_plot = dff_psth
x_plot = psth_plot.x;
h = figure;
hold on;
for i =1:size(psth_plot.mat,3)
    plot(x_plot, psth_plot.mat(:,1,i)); 
end
set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('dff'); 
title('dff locked to stim at target (stim at t=0)')
vline(0)
plot(x_plot, psth_plot.mean, 'k', 'LineWidth', 1.3)

if save_bool
    filename = 'psth_dff_traces_target_roi'
    save_ext(save_dir, filename, ext_cell)
end

%%
%--------------------------------------------------------------------------
%Target ROI PSTH + SEM

y_plot = psth_plot.mean; 
y_sem = psth_plot.sem; %-min(y_plot); 
h = figure;
hold on;

bot = y_plot-y_sem; 
bot = bot(:)'; 
top = y_plot+y_sem; 
top = top(:)'; 
inBetween = [bot fliplr(top)];
fill([x_plot, fliplr(x_plot)] , inBetween, 0.8*[1 1 1]); 
vline(0)

set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('DFF'); 
title('Target ROI DFF PSTH Mean + SEM'); 

if save_bool
    filename = 'psth_dff_sem_target_roi'
    save_ext(save_dir, filename, ext_cell)
end

%%
%--------------------------------------------------------------------------
%F Target ROI PSTH + Individual Traces
psth_plot = f_psth

h = figure;
hold on;
x_plot = psth_plot.x; 
for i =1:size(psth_plot.mat,3)
    plot(x_plot, psth_plot.mat(:,i)); 
end
plot(x_plot, psth_plot.mean, 'k', 'LineWidth', 1.3)

set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('f'); 
title('f locked to stim at target (stim at t=0)')
vline(0)

if save_bool
    filename = 'psth_f_traces_target_roi'
    save_ext(save_dir, filename, ext_cell)
end

%%
%--------------------------------------------------------------------------
%Target ROI PSTH + SEM
psth_plot = f_psth

y_plot = psth_plot.mean; 
y_sem = psth_plot.sem; %-min(y_plot); 
h = figure;
hold on;

bot = y_plot-y_sem; 
bot = bot(:)'; 
top = y_plot+y_sem; 
top = top(:)'; 
inBetween = [bot fliplr(top)];
fill([x_plot, fliplr(x_plot)] , inBetween, 0.8*[1 1 1]); 
vline(0)

set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('DFF'); 
title('Target ROI F PSTH Mean + SEM'); 

if save_bool
    filename = 'psth_f_sem_target_roi'
    save_ext(save_dir, filename, ext_cell)
end

%%
%--------------------------------------------------------------------------
%Grid PSTH
psth_win = grid_psth_win; 
x_plot = (psth_win(1):psth_win(2))/30; 

incr_width = 1/(num_x_incr+1); 
incr_height = 1/(num_y_incr+1); 

%Range can be changed: 
ymin = -0.1;
ymax = 1.5; 

plot_individual_trials = 1; 

h = figure;%
h = figure('Position', [screensize(3)/2 1 screensize(3)/4 screensize(3)/4]);%divide by 2.5 to have gray space
for grid_idx = 1:num_grid_pts %num_grid_pts%1:num_grid_pts
    ax = axes('Position', [x_grid_norm(grid_idx) y_grid_norm(grid_idx) incr_width, incr_height]); 
    hold on; 
    plot([0 0], [ymin ymax], 'r'); 
    
    if plot_individual_trials
        mat_plot = response_over_grid(grid_idx).psth_mat; 
        num_trials = size(mat_plot,3); 
        for j = 1:num_trials
            plot(x_plot, mat_plot(:,1,j), 'Color', 0.8*ones(3,1)); 
        end
        plot(x_plot, response_over_grid(grid_idx).psth_mean, 'k'); 
    end    
    ylim([ymin ymax]); 
    set(gca, 'xticklabel',[])
    set(gca, 'yticklabel',[])
%     axis square
end

if save_bool
    filename = 'grid_psth'
    save_ext(save_dir, filename, ext_cell)
end

%%
%change in dff: 
%--------------------------------------------------------------------------
scatter_size = 400; 
h = figure; 
imagesc(target_roi_data.im_bg);
hold on; 
for grid_idx = 1:num_grid_pts
    val = change_over_grid(grid_idx).delta_mean;
%     scatter(x_mesh_flat(grid_idx), y_mesh_flat(grid_idx), scatter_size, val, 'filled', 'MarkerFaceAlpha', 0.3); 
    scatter(x_mesh_flat(grid_idx), y_mesh_flat(grid_idx), scatter_size, val, 'x', 'LineWidth', 1.5); 
end
colormap('jet'); 
% caxis([0 2000]); 
colorbar
axis square
set(h, 'Position', screensize); 

% title('Percent change in DFF of target ROI for stims marked on image'); 
title('change in DFF of target ROI for stims marked on image'); 
if save_bool
    filename = 'stim_response_on_bg_change'
    save_ext(save_dir, filename, ext_cell)
end


%%
%perc change in dff: 
%--------------------------------------------------------------------------
scatter_size = 400; 
h = figure; 
imagesc(target_roi_data.im_bg);
hold on; 
for grid_idx = 1:num_grid_pts
    val = change_over_grid(grid_idx).perc_change_mean;
    scatter(x_mesh_flat(grid_idx), y_mesh_flat(grid_idx), scatter_size, val, 'x', 'LineWidth', 1.5); 
end
colormap('jet'); 
caxis([0 2000]); 
colorbar
axis square
set(h, 'Position', screensize); 

title('Percent change in DFF of target ROI for stims marked on image'); 
if save_bool
    filename = 'stim_response_on_bg_perc_change'
    save_ext(save_dir, filename, ext_cell)
end

%%
%Response over Distance: 
%Make a vector of stim distance to target point, and the percent change of
%response, trial averaged
d_px    = []; 
d_um    = [];
r_delta = []; 
r_perc  = [];
%Target: x = col, y = row
for grid_idx = 1:num_grid_pts
    d_px_i = sqrt((response_over_grid(grid_idx).col - target_roi.x)^2 + (response_over_grid(grid_idx).row - target_roi.y)^2);
    d_um_i = sqrt(micronsPerPixel.x^2*(response_over_grid(grid_idx).col - target_roi.x)^2 + micronsPerPixel.y^2*(response_over_grid(grid_idx).row - target_roi.y)^2);
    d_px = [d_px d_px_i]; 
    d_um = [d_um d_um_i]; 
    
    r_delta_i   = change_over_grid(grid_idx).delta_mean;
    r_delta     = [r_delta r_delta_i]; 
    
    r_perc_i    = change_over_grid(grid_idx).perc_change_mean;
    r_perc      = [r_perc r_perc_i]; 
end

%%
% h = figure;
% scatter(d_px,r_delta, 'fill'); 
% % xlim([-10 170]); 
% xlim([-10 80]);
% xlabel('distance between grid stim point and target roi (pixels)'); 
% ylabel('change in DFF, 1s before and after stim'); 
% axis square
% set(gca,'TickDir','out');
% if save_bool
%     filename = 'response_change_vs_distance_px'
%     save_ext(save_dir, filename, ext_cell)
% end
% 
% h = figure;
% scatter(d_px,r_perc, 'fill'); 
% % xlim([-10 170]); 
% xlim([-10 80]);
% ylim([-100 100]); 
% xlabel('distance between grid stim point and target roi (pixels)'); 
% ylabel('perc change in DFF, 1s before and after stim'); 
% axis square
% set(gca,'TickDir','out');
% if save_bool
%     filename = 'response_perc_vs_distance_px'
%     save_ext(save_dir, filename, ext_cell)
% end

%
% h = figure;
% scatter(d_um,r_delta, 'fill'); 
% xlim([0 60]); 
% xlabel('distance between grid stim point and target roi (um)'); 
% ylabel('change in DFF, 1s before and after stim'); 
% axis square
% set(gca,'TickDir','out');
% if save_bool
%     filename = 'response_change_vs_distance_um'
%     save_ext(save_dir, filename, ext_cell)
% end
% 
% h = figure;
% scatter(d_um,r_perc, 'fill'); 
% % xlim([-10 170]); 
% xlim([0 60]); 
% ylim([-100 2000]); 
% xlabel('distance between grid stim point and target roi (um)'); 
% ylabel('change in DFF, 1s before and after stim'); 
% axis square
% set(gca,'TickDir','out');
% if save_bool
%     filename = 'response_perc_vs_distance_um'
%     save_ext(save_dir, filename, ext_cell)
% end

% % save_path = fullfile(save_dir, 'plot_data.mat'); 
% % save(save_path); 
end

function save_ext(save_dir, filename, ext_cell)
    for i = 1:length(ext_cell)
        ext_i = ext_cell{i};
        save_path = fullfile(save_dir, [filename ext_i]); 
        export_fig(save_path); 
    end  
end