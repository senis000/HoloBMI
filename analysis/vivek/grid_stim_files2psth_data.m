function grid_stim_files2psth_data(psth_win_samples, stim_lag_frames, tif_dir, grid_data_path, im_bg_path, v_path, frame_chan, stim_chan, save_path)
%stim_lag: lag (num frames) between frame at which stim occurred in voltage
%data and when it actually occurred.  This lag can occur if the imaging
%session drops data :( 
%
% In experiments from Feb-Mar 2020:
% frame_chan = 1+1+1;
% stim_chan = 1+6+1; 

stim_lag_time = stim_lag_frames/30;
load(grid_data_path) %contains data about the grid for stimulation
%x_mesh_flat, y_mesh_flat have the idxs for grid points
%ctr_idx: idx for the ctr point
num_stims = length(stim_sequence)
num_target_stims = length(find(stim_sequence == ctr_idx))
im_bg = imread(im_bg_path); 
im_bg_min_perc = 0; 
im_bg_max_perc = 99.9; 
[im_bg, im_bg_min, im_bg_max] = scale_im(im_bg, im_bg_min_perc, im_bg_max_perc);

%--------------------------------------------------------------------------
%Grid Variables
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

%--------------------------------------------------------------------------
%Mask containing all grid ROI
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
chan_data = struct(...
    'label', 'g', ...
    'chan_idx', 2); %in RGB, G is 2nd
grid_roi_data = label_mask2roi_data_single_channel(im_bg, grid_mask, chan_data);
%This produces plots...

%--------------------------------------------------------------------------
%Mask containing target ROI
target_roi.x = target.x;
target_roi.y = target.y; 
target_roi.r = 10; 
%The following always works, but doesn't use roi drawn online.  
%circle_mask(dim, target_roi.y, target_roi.x, target_roi.r); 

%The following works in later experiments:
target_mask = target_roi_data.roi_bin_cell{1}; 

%--------------------------------------------------------------------------
%Load the voltageRec file: 
disp('loading voltageRec...')
tic
voltageRec = readmatrix(v_path);
toc
disp('loaded voltageRec!'); 

%--------------------------------------------------------------------------
%Extract, normalize the frame and stim signals
min_percentile = 0; 
max_percentile = 100; 
[frame_pulse, frame_min, frame_max, frame_range] = range_norm(voltageRec(:, frame_chan), max_percentile, min_percentile);

min_percentile = 0; 
max_percentile = 100; 
[stim_pulse, stim_min, stim_max, stim_range] = range_norm(voltageRec(:, stim_chan), max_percentile, min_percentile);

%--------------------------------------------------------------------------
%Find the edges in frame, stim signals
%Map the stim_idx to the frame_idx that preceds it
edge_thresh_coeff = 0.5; 
stim_edge   = find_edges_in_pulse_data(stim_pulse, edge_thresh_coeff);
frame_edge  = find_edges_in_pulse_data(frame_pulse, edge_thresh_coeff);
%Map the stim_idx to the frame_idx that precedes 
stim_frame = []; %zeros(size(stim_edge.rise));
for i=stim_edge.rise(:)'
    i_diff = i-frame_edge.rise; 
    i_diff(i_diff > 0) = -1e7; 
    [M,I] = max(i_diff);
    stim_frame = [stim_frame I]; 
end
stim_frame = stim_frame-stim_lag_frames;
target_stim_frame = stim_frame(find(stim_sequence == ctr_idx))-1;

%--------------------------------------------------------------------------
%Extract the target ROI's F from the movie. 
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

%--------------------------------------------------------------------------
%Calculate dff: sliding window f0
%some standard parameters: 2 minute sliding window F0, F0 percentile = 10
f0_win = 2*60*30; 
f0_perc = 10; 
ts = target_f_ts.'; 
[dff, f0] = calc_dff_f0_sliding_perc(ts, f0_perc, f0_win);

%--------------------------------------------------------------------------
%Range normalize the f, dff for plotting / analysis: 
[dff_norm,~] = range_norm(dff, 100, 1);
[f_norm,~] = range_norm(ts', 100, 1); 

%--------------------------------------------------------------------------
%f psth: 
psth_win = psth_win_samples; 
% psth_win = 30*[-8 8];
ts = target_f_ts.'; 
[psth_mean, psth_sem, stim_psth_mat] = calc_psth(ts, target_stim_frame, psth_win);

f_psth.mean         = psth_mean; 
f_psth.sem          = psth_sem; 
f_psth.mat          = stim_psth_mat;
f_psth.x            = (psth_win(1):psth_win(2))/30;
f_psth.psth_win     = psth_win; 

% psth_win = 30*[-8 8];
[psth_mean, psth_sem, stim_psth_mat] = calc_psth(dff, target_stim_frame, psth_win);
dff_psth.mean       = psth_mean; 
dff_psth.sem        = psth_sem; 
dff_psth.mat        = stim_psth_mat; 
dff_psth.x          = (psth_win(1):psth_win(2))/30;
dff_psth.psth_win   = psth_win; 

%--------------------------------------------------------------------------
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

grid_psth_win = 30*[-5 5];
psth_win = grid_psth_win; 
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

%--------------------------------------------------------------------------
%Change in activity over grid
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

%--------------------------------------------------------------------------
%Response as a function of distance: 
%
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
save(save_path); 