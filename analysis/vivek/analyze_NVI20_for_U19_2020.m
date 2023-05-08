%{
code taken from 'analyze_NY28_D3'

%analyze_NVI20_for_U19_2020.m
%6/30/2020:
%Goals:
%1) stim PSTH + FOV with ROI marked

Data directory:
D:\Dropbox\Data\holobmi_nuria\NVI20_D02

Files: 

BASELINE: 
BaselineOnline191106T192348.mat

STIM pretrain: 
BMI_online191106T203409.mat

BMI: 
BMI_online191106T211732.mat

%target_calibration_ALL:
target_calibration_ALL_20191106T194306.mat

%BMI_target_info: 
BMI_target_info_20191106T194306

%}

%%
%Code Paths
%--------------------------------------------------------------------------
dropbox_dir = 'D:\Dropbox\Code'
addpath(genpath(fullfile(dropbox_dir, 'analysis_util')))
addpath('D:\Dropbox\Code\export_fig\export_fig_11.19.15'); 
addpath('D:\Dropbox\Code\holobmi_git\HoloBMI/baseline_target_calibration')
addpath('D:\Dropbox\Code\holobmi_git\HoloBMI\analysis\vivek'); 

%%
%Data Paths
%--------------------------------------------------------------------------
data_dir = 'D:\Dropbox\Data\holobmi_nuria\NVI20_D02'
exist(data_dir)

%Save Plots:
plot_dir = fullfile(data_dir, 'analysis_plots'); 
mkdir(plot_dir)

im_mean_path = fullfile('C:\DATA\holobmi\NVI20_D02\baseline', 'AVG_baseline_moco.tif'); 
mask_path = fullfile(data_dir, 'strcMask.mat'); 

im_file_cell = {...
    fullfile(data_dir, 'roi_data.mat')};
% roi_data.mat
% strcMask.mat

base_file_cell = {fullfile(data_dir, 'BaselineOnline191106T192348.mat')}; 
target_info_cell = {...
    fullfile(data_dir, 'BMI_target_info_20191106T194306.mat')}; 
for i =1:length(target_info_cell)
    exist(target_info_cell{i})
end

cal_cell = {...
    fullfile(data_dir, 'target_calibration_ALL_20191106T194306.mat')}; 
for i =1:length(cal_cell)
    exist(cal_cell{i})
end

bmi_file_cell = {...
    fullfile(data_dir, 'BMI_online191106T211732.mat')}; 
for i =1:length(bmi_file_cell)
    exist(bmi_file_cell{i})
end

pretrain_file_cell = {...
    fullfile(data_dir, 'BMI_online191106T203409.mat')}; 
for i=1:length(pretrain_file_cell)
    exist(pretrain_file_cell{i})
end



%%
%Load data from files
%--------------------------------------------------------------------------
%im, base_raw, tinfo, base, bmi, pretrain
% - each is a struct array, dimension: number of files
clc
%im
for i=1:length(im_file_cell)
    if(~exist(im_file_cell{i}))
        disp(['missing im file ' num2str(i)]); 
    else
        im(i) = load(im_file_cell{i}); 
    end
end
%Image data: 
im_load = imread(im_mean_path); 
mask_data = load(mask_path); 

%base_raw
for i=1:length(base_file_cell)
    if(~exist(base_file_cell{i}))
        disp(['missing base file ' num2str(i)]); 
    else
        base_raw(i) = load(base_file_cell{i}); 
    end
end

%tinfo
for i=1:length(target_info_cell)
    if(~exist(target_info_cell{i}))
        disp(['missing target info file ' num2str(i)]);         
    else
        tinfo(i) = load(target_info_cell{i}); 
    end
end

%base - calibration
for i=1:length(cal_cell)
    if(~exist(cal_cell{i}))
        disp(['missing calibraion file ' num2str(i)]);         
    else
        base(i) = load(cal_cell{i}); 
    end
end

%bmi
for i=1:length(bmi_file_cell)
    if(~exist(bmi_file_cell{i}))
        disp(['missing bmi file ' num2str(i)]);         
    else
        bmi(i) = load(bmi_file_cell{i}); 
    end
end

%pretrain
for i=1:length(pretrain_file_cell)
    if(~exist(pretrain_file_cell{i}))
        disp(['missing pretrain file ' num2str(i)]);    
    else
        pretrain(i) = load(pretrain_file_cell{i}); 
    end
end

%%
%Load background_image:

im_double = double(im_load);
im_n = im_double/max(max(im_double)); 
min_perc = 0.1; 
max_perc = 99.7; 
[im_s, min_val_return, max_val_return] = scale_im(im_n, min_perc, max_perc);
h = figure;
imshow(im_s); 
% export_fig(h, fullfile(plot_dir, ['fov' '.eps'])); 
% export_fig(h, fullfile(plot_dir, ['fov' '.png'])); 

%%
%Show the BMI cells: 
h = figure;
imshow(im.roi_data.im_roi); 

%%
mask_data
mask_data.E_id

%%
h = figure;
imshow(im.roi_data.roi_mask); 
colormap('parula'); 
caxis([0 31])
colorbar

%%
%Process pretrain, bmi (base already processed through calibration code)
%--------------------------------------------------------------------------

%Parameters from baseline calibration: 
%baseFrames = 2*60*30; %(2 minutes with our frame rate)
%movingAverageFrames = 4; (from this version of BMI code)
%ToDo: make these more explicitly passed from baseline to the BMI code.

%---------------------------------------
% sec_per_reward_range = [120 100]; 
% 
% frames_per_reward_range = sec_per_reward_range*baseline_frameRate %[1 1.5]*60*frameRate
% %multiply by frames per minute to convert
% %to 
% 
% frameRate = 29.989;
% target_on_cov_bool = 0
% prefix_win = 40
% f0_win_bool = 1
% f0_win = 2*60*ceil(frameRate)
% dff_win_bool = 1
% dff_win = 4
%  
% reward_per_frame_range = 1./frames_per_reward_range
% 
% cursor_zscore_bool = 0;
% f0_init_slide = 0; 
% 
%  baseline2target_vBMI(n_f_file, A_file, onacid_bool,  ...
%     E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
%     prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
%     cursor_zscore_bool, f0_init_slide);
%---------------------------------------

movingAverageFrames = 4; 
cursor_zscore_bool = 0; 

%%
%base_est, pretrain_est, bmi_est
%
%collect data into format for analysis.  
%reconstruct cursor from neural activity for confirmation

%Baseline:
%Push baseline data through same function, for parallel structure for plotting
clear base_est
disp('processing baseline'); 
for i = 1:length(base)
    f_preprocessed  = 1; 
    f0_i            = base(i).f0.';
    f_i             = base(i).f_smooth';
    bData_i         = tinfo(i); %bData for baseline
    cursor_i        = base(i).cursor_obs.'; 
    base_est(i) = ...
        base_bmiAct2cursor(f_preprocessed, f0_i, f_i, movingAverageFrames, cursor_zscore_bool, bData_i, cursor_i);
    disp('number of errors in cursor recon:'); 
    sum(base_est(i).cursor - base_est(i).data_valid.cursor > 1e-3)  
    
    %Pads data to original size: 
%     base_est(i) = ...
%         base_bmiAct2cursor_v2(f_preprocessed, f0_i, f_i, movingAverageFrames, cursor_zscore_bool, bData_i, cursor_i);    
%     
%     disp('number of errors in cursor recon:'); 
%     valid = base_est(i).c_valid_idxs; 
%     sum(base_est(i).cursor(valid) - base_est(i).data_valid.cursor > 1e-3)    
end

%%
%Pretrain:
clear pretrain_est
disp('processing pretrain'); 
for i=1:length(pretrain)
    f_preprocessed  = 0; 
    f0_i        = pretrain(i).data.baseVector;
    f_i         = pretrain(i).data.bmiAct;
    bData_i     = pretrain(i).bData; 
    cursor_i    = pretrain(i).data.cursor;   
    
    pretrain_est(i) = ...
        base_bmiAct2cursor(f_preprocessed, f0_i, f_i, movingAverageFrames, cursor_zscore_bool, bData_i, cursor_i);
    disp('number of errors in cursor recon:'); 
    sum(pretrain_est(i).cursor - pretrain_est(i).data_valid.cursor > 1e-3)
    
%     pretrain_est(i) = ...
%         base_bmiAct2cursor_v2(f_preprocessed, f0_i, f_i, movingAverageFrames, cursor_zscore_bool, bData_i, cursor_i);
%     disp('number of errors in cursor recon:'); 
%     valid = pretrain_est(i).c_valid_idxs; 
%     sum(pretrain_est(i).cursor(valid) - pretrain_est(i).data_valid.cursor > 1e-3)
end

%%

%BMI: 
clear bmi_est
disp('processing bmi'); 
for i=1:length(bmi)
    f_preprocessed  = 0; 
    f0_i        = bmi(i).data.baseVector;
    f_i         = bmi(i).data.bmiAct;
    bData_i     = bmi(i).bData; 
    cursor_i    = bmi(i).data.cursor;    

    bmi_est(i) = ...
        base_bmiAct2cursor(f_preprocessed, f0_i, f_i, movingAverageFrames, cursor_zscore_bool, bData_i, cursor_i);
    disp('number of errors in cursor recon:'); 
    sum(bmi_est(i).cursor - bmi_est(i).data_valid.cursor > 1e-3)
    
%     bmi_est(i) = ...
%         base_bmiAct2cursor_v2(f_preprocessed, f0_i, f_i, movingAverageFrames, cursor_zscore_bool, bData_i, cursor_i);
%     disp('number of errors in cursor recon:');
%     valid = bmi_est(i).c_valid_idxs; 
%     sum(bmi_est(i).cursor(valid) - bmi_est(i).data_valid.cursor > 1e-3)    
end

%%
%Plot co-expression and BMI E1, E2 cells
%--------------------------------------------------------------------------
%ToDo formally: 



%%
%Plot each block's activity.  Plot f and f0 overlaid, plot dff
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
E_id = [ones(1,4) 2*ones(1,4)]'; 
num_BMI = length(tinfo.E1_base) + length(tinfo.E2_base);

%%
%BASE: 
%--------------------------------------------------------------------------
frame_rate  = 30; 
time_scale  = 1/(frame_rate*60); 
block_str   = 'base'; 
block_data  = base_est; 

%Define Events to Lock to: 
clear event_struct
for i = 1:length(block_data)
    valid_idxs                  = block_data(i).data_valid.valid_idxs;
    
    event_i = 0; 
    %Target Hit: 
    event_i     = event_i +1;
    hit_idxs    = base(i).valid_hit_idxs; 
    hit_time                    = time_scale*hit_idxs; 
    
    event_struct(i,event_i).data      = hit_time;     
    event_struct(i,event_i).label     = 'hit'; %hit_b2b'; 
    event_struct(i,event_i).valid     = 1;     
    event_struct(i,event_i).line      = '--k';
end
base_event = event_struct; 
%%
win_zoom    = [0 5]; 
save_dir    = plot_dir; 
save_bool   = 0;  
plot_big    = 0; 
plot_block_neural(block_str, block_data, event_struct, frame_rate, win_zoom, save_dir, save_bool, plot_big)


%%
%PRETRAIN: 
%--------------------------------------------------------------------------
frame_rate  = 30; 
time_scale  = 1/(frame_rate*60); 
block_str   = 'pretrain'; 
block_data  = pretrain_est; 

%Define Events to Lock to: 
for i = 1:length(pretrain_est)
    valid_idxs                  = pretrain_est(i).data_valid.valid_idxs;
    
    event_i = 0; 
    %Target Hit: 
    event_i     = event_i +1;
    
    hit_idxs    = find(pretrain(i).data.selfHits(valid_idxs)); 
    cursor      = pretrain_est(i).data_valid.cursor; 
    base_val    = pretrain_est(i).bData.T1/2; 
    b2base_num_samples = 1; 
    
    length(hit_idxs)
    [valid_hit_idxs] = ...
        analyze_hits_from_cursor_base(cursor, hit_idxs, base_val, b2base_num_samples);    
    length(valid_hit_idxs)
    
    hit_time                    = time_scale*valid_hit_idxs; 
    
    event_struct(i,event_i).data      = hit_time;     
    event_struct(i,event_i).label     = 'hit_b2b'; 
    event_struct(i,event_i).valid     = 1;     
    event_struct(i,event_i).line      = '--k';
    
    %STIM:
    event_i     = event_i +1;    
    
    valid_idxs                  = pretrain_est(i).data_valid.valid_idxs;
    stim_time                   = time_scale*find(pretrain(i).data.holoDelivery(valid_idxs)); 
    event_struct(i,event_i).data      = stim_time;     
    event_struct(i,event_i).label     = 'stim'; 
    event_struct(i,event_i).valid     = 1; 
    event_struct(i,event_i).line      = '--r'; 
    
    %STIM HIT:
    event_i     = event_i +1;
    
    valid_idxs                  = pretrain_est(i).data_valid.valid_idxs;
    stim_time                   = time_scale*find(pretrain(i).data.holoHits(valid_idxs)); 
    event_struct(i,event_i).data      = stim_time;     
    event_struct(i,event_i).label     = 'stim_hit'; 
    event_struct(i,event_i).valid     = 1;
    event_struct(i,event_i).line      = '--b'; 
    
end
pretrain_event = event_struct; 

%%
%TODO: fix given the way i changed preprocessing... NaN is messing up the
%plot

close all; 
win_zoom    = [0 5]; 
save_dir    = plot_dir; 
save_bool   = 0; 
plot_big    = 0; 
plot_block_neural(block_str, block_data, event_struct, frame_rate, win_zoom, save_dir, save_bool, plot_big)


%%
%Collect BMI neural timeseries: 
%num_neurons x num_samples

clear n
n.base          = base.n_analyze.'; %which should be dff
n.pretrain      = pretrain_est.dff; 
n.bmi           = bmi_est.dff;

n.base_z        = base.dff_z.'; 
n.pretrain_z    = pretrain_est.dff_z; 
n.bmi_z         = bmi_est.dff_z;

n_pad.base      = base_est.pad.dff;
n_pad.pretrain  = pretrain_est.pad.dff; 
n_pad.bmi       = bmi_est.pad.dff; 

n_pad.base_z      = base_est.pad.dff_z;
n_pad.pretrain_z  = pretrain_est.pad.dff_z; 
n_pad.bmi_z       = bmi_est.pad.dff_z; 

%%
% Collect cursor:
clear cursor
%base:
cursor.base = base.cursor_obs.'; 
%pretrain:
valid_idxs = find(~isnan(pretrain.data.cursor));
cursor.pretrain = pretrain.data.cursor(valid_idxs); 
%bmi: 
valid_idxs = find(~isnan(bmi.data.cursor));
cursor.bmi = bmi.data.cursor(valid_idxs); 
cursor

min_data = min([cursor.base(:); cursor.pretrain(:); cursor.bmi(:)]); 
max_data = max([cursor.base(:); cursor.pretrain(:); cursor.bmi(:)]); 
num_bins = 51; 
num_bin_edges = num_bins-1; 
bin_edges = linspace(min_data, max_data, num_bin_edges); 
[bin_centers] = bin_edges2bin_centers(bin_edges);

%%
[data_binned, base_bin] = bin_data(bin_edges, cursor.base);
base_bin_n = base_bin/sum(base_bin); 

[data_binned, pretrain_bin] = bin_data(bin_edges, cursor.pretrain);
pretrain_bin_n = pretrain_bin/sum(pretrain_bin); 

[data_binned, bmi_bin] = bin_data(bin_edges, cursor.bmi);
bmi_bin_n = bmi_bin/sum(bmi_bin); 

% plot(x_data, y_data, '.-', 'LineWidth', 1.5, 'MarkerSize', 10, 'color', colors(session_idx+color_offset, :))

%%
%Unnormalized histogram:
h = figure;
hold on; 
plot(bin_centers, base_bin, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, pretrain_bin, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, bmi_bin, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
legend({'base', 'pretrain', 'bmi'}); 
xlim([-1.5 1]); 
title('baseline'); 
%TODO: get hist code I used when doing fernandoBMI analysis

%%
%Normalized histogram
h = figure;
hold on; 
plot(bin_centers, base_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, pretrain_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, bmi_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
vline(base.T, 'r:'); 
legend({'base', 'pretrain', 'bmi', 'target'}); 

xlim([-1.5 1]); 
xlabel('cursor value: E2-E1'); 
ylabel('frac time entered'); 
set(gca,'TickDir','out');
% export_fig(h, fullfile(plot_dir, 'base_pretrain_bmi_hist.eps')); 
%TODO: get hist code I used when doing fernandoBMI analysis
%%
h = figure;
hold on; 
plot(bin_centers, base_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, pretrain_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, bmi_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
vline(base.T, 'r:'); 
legend({'base', 'pretrain', 'bmi', 'target'}); 

xlim([0 1]); 
xlabel('cursor value: E2-E1'); 
ylabel('frac time entered'); 
set(gca,'TickDir','out');
% export_fig(h, fullfile(plot_dir, 'base_pretrain_bmi_hist_zoom.eps')); 
%TODO: get hist code I used when doing fernandoBMI analysis

%%
%Take cursor data after stim: 
%Hist: 
holo_valid = find(pretrain.data.holoDelivery);
holo_win = [];
num_frames_post_stim = 30; 
for i=1:length(holo_valid)
    holo_win = [holo_win holo_valid(i):(holo_valid(i)+num_frames_post_stim)]; 
end
holo_cursor = pretrain.data.cursor(holo_win); 
% holo_cursor = pretrain_est.pad.cursor(holo_win);
[data_binned, holo_cursor_bin] = bin_data(bin_edges, holo_cursor);
holo_cursor_bin_n = holo_cursor_bin/sum(holo_cursor_bin); 

%%
%data outside of stim 

%%
h = figure;
hold on; 
plot(bin_centers, pretrain_bin_n, '.-k', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, holo_cursor_bin_n, '.-r', 'LineWidth', 1.5, 'MarkerSize', 10); 
xlim([-2 2]); 
vline(base.T)
legend({'pretrain E2-E1', 'stim E2-E1'}); 
xlabel('E2-E1'); 
ylabel('frac time entered'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_cursor.eps')); 


%%
n.pretrain_holo = n_pad.pretrain(:, holo_win); 

%%
n.base_E1 = sum(n.base(1:4,:), 1); 
n.base_E2 = sum(n.base(5:8,:), 1); 

n.pretrain_E1 = sum(n.pretrain(1:4,:), 1); 
n.pretrain_E2 = sum(n.pretrain(5:8,:), 1); 

%
n.pretrain_holo_E1 = sum(n.pretrain_holo(1:4,:), 1); 
n.pretrain_holo_E2 = sum(n.pretrain_holo(5:8,:), 1); 

%
%E1: 
min_data = min([n.pretrain_E1 n.pretrain_holo_E1]); 
max_data = max([n.pretrain_E1 n.pretrain_holo_E1]); 

% min_data = min([n.base_E1 n.pretrain_holo_E1]); 
% max_data = max([n.base_E1 n.pretrain_holo_E1]); 

num_bins = 25; 
num_bin_edges_E1 = num_bins-1; 
bin_edges_E1 = linspace(min_data, max_data, num_bin_edges_E1); 
[bin_centers_E1] = bin_edges2bin_centers(bin_edges_E1);

[data_binned, base_E1_bin] = bin_data(bin_edges_E1, n.base_E1);
[data_binned, pretrain_E1_bin] = bin_data(bin_edges_E1, n.pretrain_E1);
[data_binned, holo_E1_bin] = bin_data(bin_edges_E1, n.pretrain_holo_E1);

base_E1_bin_n = base_E1_bin/sum(base_E1_bin); 
pretrain_E1_bin_n = pretrain_E1_bin/sum(pretrain_E1_bin); 
holo_E1_bin_n = holo_E1_bin/sum(holo_E1_bin); 


%%

h = figure;
hold on; 
plot(bin_centers_E1, pretrain_E1_bin_n, '.-k', 'LineWidth', 1.5, 'MarkerSize', 10); 
% plot(bin_centers_E1, base_E1_bin_n, '.-k', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers_E1, holo_E1_bin_n, '.-r', 'LineWidth', 1.5, 'MarkerSize', 10); 
% xlim([-1 1]); 
% legend({'baseline E1', 'stim E1'}); 
legend({'pretrain E1', 'stim E1'}); 
xlabel('E1'); 
ylabel('frac time entered'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_E1.eps')); 
export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_E1.png')); 

%%
%E2: 
% min_data = min([n.base_E2 n.pretrain_holo_E2]); 
% max_data = max([n.base_E2 n.pretrain_holo_E2]); 
min_data = min([n.pretrain_E2 n.pretrain_holo_E2]); 
max_data = max([n.pretrain_E2 n.pretrain_holo_E2]); 

num_bins = 30; 
num_bin_edges_E2 = num_bins-1; 
bin_edges_E2 = linspace(min_data, max_data, num_bin_edges_E2); 
[bin_centers_E2] = bin_edges2bin_centers(bin_edges_E2);


% [data_binned, base_E2_bin] = bin_data(bin_edges_E2, n.base_E2);
[data_binned, pretrain_E2_bin] = bin_data(bin_edges_E2, n.pretrain_E2);
[data_binned, holo_E2_bin] = bin_data(bin_edges_E2, n.pretrain_holo_E2);

% base_E2_bin_n = base_E2_bin/sum(base_E2_bin); 
pretrain_E2_bin_n = pretrain_E2_bin/sum(pretrain_E2_bin); 
holo_E2_bin_n = holo_E2_bin/sum(holo_E2_bin); 

%%
h = figure;
hold on; 
% plot(bin_centers_E2, base_E2_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers_E2, pretrain_E2_bin_n, '.-k', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers_E2, holo_E2_bin_n, '.-r', 'LineWidth', 1.5, 'MarkerSize', 10); 
% xlim([-1 1]); 
legend({'pretrain E2', 'stim E2'}); 
xlabel('E2'); 
ylabel('frac time entered'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_E2.eps')); 
export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_E2.png')); 

%%
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 

% n_plot = zscore(n.base.'); 
n_plot = n.base.'; 
t = (1:size(n_plot,1))/(26*60); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color);
title('base'); 

set(gca,'TickDir','out');
% export_fig(h, fullfile(plot_dir, 'base.eps')); 


%
n_plot = n.pretrain.'; 
t = (1:size(n_plot,1))/(30*60); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color);
title('pretrain'); 

%
n_plot = n.bmi.'; 
t = (1:size(n_plot,1))/(30*60); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color);
title('bmi'); 

%%
data_mat = n.pretrain.'; 
event_idxs = find(pretrain.data.holoDelivery(pretrain_est.data_valid.valid_idxs))+1;
win = [-150 150]; 
win_len = win(2)-win(1); 
[psth_mean, psth_sem, psth_mat] = calc_psth(data_mat, event_idxs, win);
% calc_psth = 
t_plot = (win(1):win(2))/30; 
h = figure; hold on;
offset = 0; 
for i=1:num_neurons
    y_plot = psth_mean(:,i); 
    y_plot = y_plot-min(y_plot);
    y_amp = max(y_plot+psth_sem(:,i)); 
    offset = offset + y_amp + 0.15; 
    y_sem = psth_sem(:,i)-min(y_plot); 
    
    plot(t_plot, y_plot-offset, 'Color', E_color{(E_id(i))}); 
    errbar(t_plot, y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
end
vline(0)
% vline((psth_win(2)-psth_win(1))/2+1); 
xlabel('time (sec)')
ylabel('deltaF over F'); 
title('PSTH of Ensemble Activity Locked to 2p stimulation'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, ['pretrain_stim_psth_win' num2str(win_len) '.eps'])); 
export_fig(h, fullfile(plot_dir, ['pretrain_stim_psth_win' num2str(win_len) '.png'])); 

%%
%Zscore, locked to stim
data_mat = zscore(n.pretrain.'); 
event_idxs = find(pretrain.data.holoDelivery(pretrain_est.data_valid.valid_idxs));
win = [-150 150]; 
win_len = win(2)-win(1); 
[psth_mean, psth_sem, psth_mat] = calc_psth(data_mat, event_idxs, win);
% calc_psth = 
%
t_plot = (win(1):win(2))/30; 
h = figure; hold on;
offset = 0; 
for i=1:num_neurons
    y_plot = psth_mean(:,i); 
    y_plot = y_plot-min(y_plot);
    y_amp = max(y_plot+psth_sem(:,i)); 
    offset = offset + y_amp + 0.15; 
    y_sem = psth_sem(:,i)-min(y_plot); 
    
    plot(t_plot, y_plot-offset, 'Color', E_color{(E_id(i))}); 
    errbar(t_plot, y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
end
vline(0)
% vline((psth_win(2)-psth_win(1))/2+1); 
xlabel('time (sec)')
ylabel('z-score deltaF over F'); 
title('PSTH of Ensemble Activity Locked to 2p stimulation'); 
set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, ['pretrain_stim_psth_Z_win' num2str(win_len) '.eps'])); 
export_fig(h, fullfile(plot_dir, ['pretrain_stim_psth_Z_win' num2str(win_len) '.png'])); 