
%%
%code directory: 
code_home = 'D:\Dropbox\Code\holobmi_git\HoloBMI'
cd(code_home)
addpath(genpath(cd))
%%
%file to load: 
clear all
date_str    = '200311'
file_analyze = 'fine_gridstim-004' %'stim-005'
%TODO: select the files to load
%fine_gridstim-000
%fine_gridstim-001
%fine_gridstim-002
%fine_gridstim-003
%fine_gridstim-004
%
%coarse_gridstim-000
%coarse_gridstim-001
%coarse_gridstim-002
%coarse_gridstim-003
%coarse_gridstim-004

%tif_dir     = ['D:\DATA\grid_stim\' date_str '\NVI20\D0\im\' file_analyze] %coarse_grid-000'
tif_dir     = 'D:\DATA\grid_stim\200311\NVI20\D0\im\fine_gridstim-004'
data_dir    = 'D:\DATA\grid_stim\200311\NVI20\D0'
load_path   = fullfile(data_dir, 'fine_grid.mat')
im_bg_path          = fullfile('D:\DATA\grid_stim\200311\NVI20\D0\redgreen', 'green.tif'); 
voltageRec_path     = ...
    fullfile('D:\DATA\grid_stim\200311\NVI20\D0\im\fine_gridstim-004', 'fine_gridstim-004_Cycle00001_VoltageRecording_001.csv');
%fullfile('D:\DATA\grid_stim\200311\NVI20\D0\im\fine_gridstim-000', 'fine_gridstim-000_Cycle00001_VoltageRecording_001.csv');  
%fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']); 

%['D:\DATA\grid_stim\' date_str '\NVI20\D0\coarse_grid_data' %['D:\DATA\grid_stim\200214\NVI20\D0\data_' file_analyze];

save_dir    = fullfile('D:\Dropbox\Data\grid_stim', date_str, file_analyze); 
mkdir(save_dir); 


%Code paths to add: 
addpath('D:\Dropbox\Code\export_fig\export_fig_3.3.20'); 

%
screensize = get(0, 'ScreenSize')

%clear
load(load_path)

%%
%This function depends on imagej taking average video

exist(im_bg_path)
im_bg = imread(im_bg_path); 

%%
h = figure;
imshow(im_bg)

im_sc_struct = struct(...
    'im', [], ...
    'minmax_perc', [], ...
    'minmax', [], ...
    'min', [], ...
    'min_perc', [], ...
    'max', [], ...
    'max_perc', []); 
num_im_sc = 0; 
[im_sc_struct, num_im_sc] = scale_im_interactive(im_bg, im_sc_struct, num_im_sc);
im_bg = im_sc_struct(end).im; 

%%
%Scatter the grid points: 

% h = figure; 
% hold on;
% scatter(x_mesh_flat, y_mesh_flat, 'x'); 
% scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 'r', 'x', 'LineWidth', 2); 
% axis square
% 
% h = figure; 
% imagesc(im_bg); 
% colormap('gray')
% hold on
% scatter(x_mesh_flat, y_mesh_flat, 'x'); 
% scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 'r', 'x', 'LineWidth', 1); 
% xlim([1 dim])
% ylim([1 dim])
% axis square

h = figure('Position', [screensize(3)/2 1 screensize(3)/2 screensize(4)]);
imagesc(im_bg); 
colormap('gray')
hold on
scatter(x_mesh_flat, y_mesh_flat, 200, 'x', 'LineWidth', 1.5); 
scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 200, 'r', 'x', 'LineWidth', 1.5); 
xlim([1 dim])
ylim([1 dim])
axis square
%Reminder, mapping from x,y to image space: 
%x = column, y = row

save_bool = 0
if save_bool
    tic
    disp('saving...')
    save_path = fullfile(save_dir, 'grid_stim.eps');
%     saveas(h, save_path); 
    export_fig(save_path); 
    
    save_path = fullfile(save_dir, 'grid_stim.png');
%     saveas(h, save_path); 
    export_fig(save_path);
    disp('DONE saving'); 
    toc
else
    disp('not saving...'); 
end
% close all; 

%%
%ROI_mask: 
grid_mask = zeros(dim, dim); 
r = 10; %To Do: figure out how big the stim radius is.  
for grid_idx = 1:num_grid_pts
    grid_col = x_mesh_flat(grid_idx); 
    grid_row = y_mesh_flat(grid_idx); 
    grid_mask_i = circle_mask(dim, grid_row, grid_col, r); 
    %This could be implemented way more efficiently.... by applying shift
    %operator in fourier domain on the circle iamge...
    grid_mask(grid_mask_i > 0) = grid_idx; 
%     [mask] =  circle_mask(dim, ctr_row, ctr_col, r)
end

debug_bool = 0
if debug_bool
    h = figure; 
    imagesc(grid_mask); 
    axis equal 
end

chan_data = struct(...
    'label', 'g', ...
    'chan_idx', 2); %in RGB, G is 2nd
grid_roi_data = label_mask2roi_data_single_channel(im_bg, grid_mask, chan_data);

%%
%Make the ROI: 
%'target' has the center, and radius
% in future we can modify script so that I draw the ROI around the target
% cell, and save it online during the grid stim expt.  (This also helps me
% do manual motion correction)

%Set manually from inspecting the video in ImageJ
% target_roi = target; %If the target data matches the neuron roi
target_roi.x = target.x;
target_roi.y = target.y; 
target_roi.r = 10; 

% target_roi.x = 252;
% target_roi.y = 280; 
% target_roi.r = 10; 
%In future we want to pick the cell, move the field of view to center,
%check the stim, then draw the roi there, then use our online monitoring of
%the cell to get the result right away.  

target_mask = target_roi_data.roi_bin_cell{1};%circle_mask(dim, target_roi.y, target_roi.x, target_roi.r); 
%%x: 248, y: 280
%
h = figure;
imagesc(target_mask)
axis square

%
% chan_data = struct(...
%     'label', 'g', ...
%     'chan_idx', 2); %in RGB, G is 2nd
% target_roi_data = label_mask2roi_data_single_channel(im_bg, target_mask, chan_data);

%%
%Number and Order of Grid Stims
save_bool = 1

num_stims = length(stim_sequence)
num_target_stims = length(find(stim_sequence == ctr_idx))
%Number of stimulations per sites:
h= figure;
hist(stim_sequence, 1:num_grid_pts);
set(gca,'TickDir','out');
xlabel('idx of grid point'); 
ylabel('number of stims to grid point'); 
title(['num stims per grid point, idx of grid pt at target roi: ' num2str(ctr_idx)])
disp('center point:')
ctr_idx 
if save_bool
    disp('saving stim_hist...');
    save_path = fullfile(save_dir, 'stim_hist.eps'); 
    export_fig(save_path); 
    
    save_path = fullfile(save_dir, 'stim_hist.png'); 
    export_fig(save_path); 
    disp('DONE saved!'); 
end

%Order of stimulation: 
h = figure; hold on; 
plot(stim_sequence, 'LineWidth', 0.5)
scatter(1:num_stims, stim_sequence)
set(gca,'TickDir','out');
xlabel('stim idx'); 
ylabel('grid idx stimmed'); 
title('order of stimulations')
if save_bool
    disp('saving stim order...'); 
    save_path = fullfile(save_dir, 'stim_seq.eps'); 
    export_fig(save_path); 
    
    save_path = fullfile(save_dir, 'stim_seq.png'); 
    export_fig(save_path); 
    disp('DONE saving')
end

%%
%Need to identify times of stims, confirm the right number occurred. 
%Load the voltageRec file:
%AI1: Frame In
%AI6: Monaco trigger rec
%AI7: Frame Trigger

%%
%Load the voltageRec file: 
disp('loading voltageRec...')
tic
%fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']); 
% exist(voltageRec_path)
voltageRec = readmatrix(voltageRec_path);
toc
disp('loaded voltageRec!'); 
%%
frame_chan = 1+1+1;
stim_chan = 1+6+1; 

plot_bool = 1
if plot_bool
    h = figure;
    plot(voltageRec(:, frame_chan)); 
    title('frame pulse'); 

    h = figure;
    plot(voltageRec(:,stim_chan))
    title('stim pulse'); 
end

min_percentile = 0; 
max_percentile = 100; 
[frame_pulse, frame_min, frame_max, frame_range] = range_norm(voltageRec(:, frame_chan), max_percentile, min_percentile);

min_percentile = 0; 
max_percentile = 100; 
[stim_pulse, stim_min, stim_max, stim_range] = range_norm(voltageRec(:, stim_chan), max_percentile, min_percentile);


%%
%Visualize stim pulse and frame pulse
plot_bool = 1;
if plot_bool
    plot_idxs = (1:2000) + 7000;
    y1 = frame_pulse(plot_idxs); 
    y2 = stim_pulse(plot_idxs);
    h = figure;
    hold on;
    plot(y1)
    plot(y2)
    legend('frame', 'stim'); 
    ylim([-0.1 1.1]); 
end
%%
edge_thresh_coeff = 0.5; 
stim_edge   = find_edges_in_pulse_data(stim_pulse, edge_thresh_coeff);
frame_edge  = find_edges_in_pulse_data(frame_pulse, edge_thresh_coeff);
% [edges] = find_edges_in_pulse_data(data, edge_thresh_coeff)

%
%Map the stim_idx to the frame_idx that precedes 
stim_frame = []; %zeros(size(stim_edge.rise));
for i=stim_edge.rise(:)'
    i_diff = i-frame_edge.rise; 
    i_diff(i_diff > 0) = -1e7; 
    [M,I] = max(i_diff);
    stim_frame = [stim_frame I]; 
end
%I visually verified that the target cell responds to the stimulation in
%imagej

%%
%Sanity: verify the timing of the pulses: 
%stim_edge.rise(1)
%stim_edge.fall(1)
%frame_edge.rise(2)
%frame_edge.fall(2)
% h =figure;
% plot(7900:8100, stim_pulse(7900:8100), '.-', 'MarkerSize', 5)
% h = figure;
% plot(20:50, frame_pulse(20:50))
%%
%how many frame pulses were counted and how many frames were in movie? 
%could this be the source of the delay? 

%%
%Load the movie: 
fileList = dir(fullfile(tif_dir, '*.ome.tif')); 
num_frames = length(fileList); 
target_f_ts = zeros(1, num_frames); 
% video_array = zeros(dim, dim, length(fileList)); 
disp('extracting target roi fluorescence...'); 
tic
for i = 1:num_frames
    path_i      = fullfile(fileList(i).folder, fileList(i).name); 
    frame_i     = double(imread(path_i)); 
    target_f_ts(i) = frame_i(:)'*target_mask(:); 
    i
end
disp('DONE!'); 
toc

%%
%calculate dff: sliding window f0
%some standard parameters: 2 minute sliding window F0, F0 percentile = 10
f0_win = 2*60*30; 
f0_perc = 10; 
ts = target_f_ts.'; 
[dff, f0] = calc_dff_f0_sliding_perc(ts, f0_perc, f0_win);

%%
[dff_norm,~] = range_norm(dff, 100, 1);
[f_norm,~] = range_norm(ts', 100, 1); 

h = figure;
hold on; 
plot(f_norm); 
plot(dff_norm); 
set(gca,'TickDir','out');
xlabel('frame'); 
ylabel('signal'); 
title('f and dff, black line = stim of target roi')
set(h, 'Position', screensize); 

lastGoodFrame = 12

target_stim_frame = stim_frame(find(stim_sequence == ctr_idx))-1;
vline(target_stim_frame, 'k')

legend('f norm','dff norm');


save_bool = 0
if save_bool
    plot_name = 'f_dff_target_stim'
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
%f psth: 
target_stim_frame = stim_frame(find(stim_sequence == ctr_idx)); %%
psth_win = 30*[-8 8];
[psth_mean, psth_sem, stim_psth_mat] = calc_psth(ts, target_stim_frame, psth_win);

f_psth.mean         = psth_mean; 
f_psth.sem          = psth_sem; 
f_psth.mat          = stim_psth_mat;
f_psth.x            = (psth_win(1):psth_win(2))/30;
f_psth.psth_win     = psth_win; 

target_stim_frame = stim_frame(find(stim_sequence == ctr_idx)); %%
psth_win = 30*[-8 8];
[psth_mean, psth_sem, stim_psth_mat] = calc_psth(dff, target_stim_frame, psth_win);
dff_psth.mean       = psth_mean; 
dff_psth.sem        = psth_sem; 
dff_psth.mat        = stim_psth_mat; 
dff_psth.x          = (psth_win(1):psth_win(2))/30;
dff_psth.psth_win   = psth_win; 

%%
psth_plot = dff_psth

h = figure;
hold on;

%
x_plot = psth_plot.x; 
for i =1:size(psth_plot.mat,3)
    plot(x_plot, psth_plot.mat(:,i)); 
end
plot(x_plot, psth_plot.mean, 'k', 'LineWidth', 1.3)
%

set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('dff'); 
title('dff locked to stim at target (stim at t=0)')
vline(0)
% vline(-8/30)

save_bool = 0
if save_bool
    plot_name = 'psth_dff_target_roi'
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
psth_plot = f_psth

h = figure;
hold on;

%
x_plot = psth_plot.x; 
for i =1:size(psth_plot.mat,3)
    plot(x_plot, psth_plot.mat(:,i)); 
end
plot(x_plot, psth_plot.mean, 'k', 'LineWidth', 1.3)
%

set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('f'); 
title('f locked to stim at target (stim at t=0)')
vline(0)

save_bool = 1
if save_bool
    plot_name = 'psth_f_target_roi'
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
%Plot mean + sem: 
y_plot  = dff_psth.mean; 
y_sem   = dff_psth.sem; %-min(y_plot); 
h = figure;
hold on;

bot = y_plot-y_sem; 
bot = bot(:)'; 
top = y_plot+y_sem; 
top = top(:)'; 
inBetween = [bot fliplr(top)];
fill([x_plot, fliplr(x_plot)] , inBetween, 0.8*[1 1 1]); 
% plot(x_plot, y_plot, 'k', 'LineWidth', 0.5)
vline(0)
set(gca,'TickDir','out');
xlabel('time (s)'); 
ylabel('DFF'); 
title('Target ROI PSTH Mean + SEM'); 

save_bool = 1
if save_bool
    plot_name = 'psth_sem_target_roi'
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
%show the responses
%show the grid points on the background image
%Convert grid points to xy and normalize
%x = col, y = dim - row + 1, width = x_step/dim, height = y_step/dim


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

debug_plot = 0
if debug_plot
    %Confirm the mapping from row,col to y,x
    grid_idxs = 1:9%:10
    close all

    h = figure;
    hold on;
    for grid_idx = grid_idxs
        scatter(x_grid_norm(grid_idx), y_grid_norm(grid_idx)); 
        xlim([0 1]); 
        ylim([0 1]); 
        axis square
    end

    h = figure; 
    hold on;
    for grid_idx = grid_idxs
        scatter(x_mesh_flat(grid_idx), y_mesh_flat(grid_idx)); 
        xlim([1 dim]); 
        ylim([1 dim]); 
        axis square
    end
end

%%
%Make a data matrix for the response locked to each stim point: 
stim_frame_over_grid = cell(num_grid_pts, 1); 
response_over_grid = repmat(...
    struct(...
    'row', [], ...
    'col', [], ...
    'x', [], ...
    'y', [], ...
    'psth_mean', [], ...
    'psth_sem', [], ...
    'psth_mat', []), [num_grid_pts, 1]); 


%%
psth_win = 30*[-5 5];
x_plot = (psth_win(1):psth_win(2))/30; 
for grid_idx = 1:num_grid_pts
    grid_stim_frame = stim_frame(find(stim_sequence == grid_idx));
    stim_frame_over_grid{grid_idx} = grid_stim_frame; 
    
    [psth_mean, psth_sem, psth_mat] = calc_psth(dff, grid_stim_frame, psth_win);
    
    response_over_grid(grid_idx).row = ...
        y_mesh_flat(grid_idx);
    response_over_grid(grid_idx).col = ...
        x_mesh_flat(grid_idx);
    
    response_over_grid(grid_idx).x = ...
        x_grid_norm(grid_idx);
    response_over_grid(grid_idx).y = ...
        y_grid_norm(grid_idx);
    
    response_over_grid(grid_idx).psth_mean = ...
        psth_mean;
    response_over_grid(grid_idx).psth_sem = ...
        psth_sem;
    response_over_grid(grid_idx).psth_mat = ...
        psth_mat;    
end

%%
debug_bool = 1
if debug_bool
    disp('-----------------------------------------------'); 
    disp('row (y)'); 
    disp('target:'); 
    target.y
    disp('response over grid:'); 
    response_over_grid(ctr_idx).row
    disp('y_mesh_flat:'); 
    y_mesh_flat(ctr_idx)
    disp('target roi (y) (manual posthoc)'); 
    target_roi.y

    %
    disp('-----------------------------------------------'); 
    disp('col (x)'); 
    disp('target:'); 
    target.x
    disp('response over grid:'); 
    response_over_grid(ctr_idx).col
    disp('x_mesh_flat:'); 
    x_mesh_flat(ctr_idx)
    disp('target roi (x) (manual posthoc)'); 
    target_roi.x
end

%%
change_over_grid = repmat(struct(...
    'pre_mat', [], ...
    'post_mat', [], ...
    'perc_change_mat', [], ...
    'perc_change', []), [num_grid_pts, 1]); 

pre_idxs    = find((x_plot >= -2 & x_plot < 0)); 
post_idxs   = find((x_plot >= 0 & x_plot <= 2));

%
for grid_idx = 1:num_grid_pts
    mat_i = response_over_grid(grid_idx).psth_mat; 
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

%%
h = figure;
hold on; 
for grid_idx = 1:num_grid_pts
%     val = change_over_grid(grid_idx).perc_change;
%     val = change_over_grid(grid_idx).perc_change_mean;
    val = change_over_grid(grid_idx).delta_mean;
    scatter(response_over_grid(grid_idx).x, response_over_grid(grid_idx).y, 1000, val,'filled'); 
end

colormap(jet);
% caxis([0 2000]); 
colorbar; 

% xlim([-0.1 1.1]); 
% ylim([-0.1 1.1]); 
axis square
% title('Delta between 2 sec after stim vs 2 sec before stim'); 
%%
% h = figure; 
% scatter(x_grid_norm, y_grid_norm); 
% xlim([0 1]); 
% ylim([0 1]); 
% axis square


%%
incr_width = 1/(num_x_incr+1); 
incr_height = 1/(num_y_incr+1); 

ymin = -0.1;
ymax = 1.5; 

plot_individual_trials = 1; 

h = figure;
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
%     axis square
end

save_bool = 1
if save_bool
    plot_name = 'grid_psth'
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
%Make a vector of stim distance to target point, and the percent change of
%response, trial averaged
d = []; 
r = []; 
%Target: x = col, y = row
for grid_idx = 1:num_grid_pts
    d_i = sqrt((response_over_grid(grid_idx).col - target_roi.x)^2 + (response_over_grid(grid_idx).row - target_roi.y)^2);
    d = [d d_i]; 
    
%     r_i = perc_change; 
    r = [r change_over_grid(grid_idx).delta_mean]; 
end

%%
h = figure;
scatter(d,r, 'fill'); 
xlim([-10 170]); 
xlabel('distance between grid stim point and target roi (pixels)'); 
ylabel('change in DFF, 1s before and after stim'); 
axis square
set(gca,'TickDir','out');

save_bool = 1
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
imagesc(grid_roi_data.im_roi)
axis square
title('Grid points overlay on the background image'); 
set(gca,'TickDir','out');

save_bool = 1
if save_bool
    plot_name = 'grid_pts_on_bg'
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
h= figure; 
imagesc(target_roi_data.im_roi)
axis square
title('Target ROI overlay on the background image'); 

save_bool = 1
if save_bool
    plot_name = 'target_on_bg'
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
% %Superimpose the grid and target: 
grid_and_target = grid_roi_data.im_roi; 
grid_and_target(:,:,3) = 0.5*grid_and_target(:,:,3); 
grid_and_target(:,:,1) = grid_and_target(:,:,1) + 0.5*target_roi_data.im_roi(:,:,3); 

h = figure; 
imagesc(grid_roi_data.im_bg);
title('background image'); 
axis square
h = figure; 
imagesc(grid_and_target); 
axis square
title('grid points (blue) and stim point (red)'); 

save_bool = 1
if save_bool
    plot_name = 'target_and_grid_on_bg'
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
% im_sc_struct = struct(...
%     'im', [], ...
%     'minmax_perc', [], ...
%     'minmax', [], ...
%     'min', [], ...
%     'min_perc', [], ...
%     'max', [], ...
%     'max_perc', []); 
% num_im_sc = 0; 
% [im_sc_struct, num_im_sc] = scale_im_interactive(im_bg, im_sc_struct, num_im_sc);
% im_bg_sc = im_sc_struct(end).im; 


%% 
% h = figure; 
% imagesc(grid_and_target);
% axis square

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

save_bool = 1
if save_bool
    plot_name = 'stim_response_on_bg'
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
%1) percent change as a function of distance from the target point
%2) a heat map of the percent change
%3) place the subplot of matlab at the appropriate coordinate


