%5.29.19
%
%{ 
Summary: 
Two baseline calibrations.  Use the second one.  
Baseline file: 
BaselineOnline190521T222151.mat

One pretrain file: 
BMI_online190521T232745

BMI split into two files: 
BMI_online190522T000029.mat
BMI_online190522T002004.mat


%File summary: 
%
%--------------------------------------------------------------------------
%Holo: 
% 8:35PM
% holoOnline190521T203550.mat
% 17.2 MB

%--------------------------------------------------------------------------
%Baseline: 
% 10:21PM
% BaselineOnline190521T222151.mat
% 36.9 MB

%--------------------------------------------------------------------------
%BMI_online:
%10:41pm
% BMI_online190521T224159
% 751 KB
%
%10:45pm
% BMI_online190521T224534.mat
% 327 KB
%
%11:27pm
% BMI_online190521T232745
% 8.7 MB
%
%11:37pm
% BMI_online190521T233746
% 85 KB
%
% 12:00AM
% BMI_online190522T000029.mat
%4.7MB
%
%12:20AM
%BMI_online190522T002004.mat
% 3.9MB

%--------------------------------------------------------------------------
%BMI_target_info: 
%1) 10:28pm
% BMI_target_info_20190521T222852.mat
%2) 10:42pm
% BMI_target_info_20190521T224235.mat

%--------------------------------------------------------------------------
%target_calibration_ALL:
%1) 10:27pm
% target_calibration_ALL_20190521T222709.mat
%2) 10:28pm
% target_calibration_ALL_20190521T222852.mat
% (probably changed task difficulty)
%3) 10:42pm
% target_calibration_ALL_20190521T224235.mat
%}
%%
%File used to seed the BMI: 
%BMI_online190521T232745
%(The second run of BMI also used the pretrain file to seed baseline
%value...?)

%%
%Code Paths
%--------------------------------------------------------------------------
addpath(genpath('/Users/vivekathalye/Dropbox/Code/analysis_util'))
addpath('/Users/vivekathalye/Dropbox/Code/export_fig/export_fig_11.19.15'); 
addpath('/Users/vivekathalye/Dropbox/Code/holobmi_git/HoloBMI/baseline_target_calibration')

%%
%Data Paths
%--------------------------------------------------------------------------
data_dir = '/Users/vivekathalye/Dropbox/Data/holo_may2019/NY28_190521'
exist(data_dir)

%Save Plots:
plot_dir = fullfile(data_dir, 'analysis_plots_v2'); 
mkdir(plot_dir)

im_file_cell = {...
    fullfile(data_dir, 'red.mat')};
base_file_cell = {fullfile(data_dir, 'BaselineOnline190521T222151.mat')}; 
target_info_cell = {...
    fullfile(data_dir, 'BMI_target_info_20190521T224235.mat')}; 
cal_cell = {...
    fullfile(data_dir, 'target_calibration_ALL_20190521T224235.mat')}; 
bmi_file_cell = {...
    fullfile(data_dir, 'BMI_online190522T000029.mat'), ...
    fullfile(data_dir, 'BMI_online190522T002004.mat')}; 
pretrain_file_cell = {...
    fullfile(data_dir, 'BMI_online190521T232745.mat')}; 
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
    % function [est] = base_bmiAct2cursor(baseVector, bmiAct, movingAverageFrames, cursor_zscore_bool, bData, cursor)    
%     disp('number of errors in cursor recon:'); 
    disp('number of errors in cursor recon:'); 
    sum(base_est(i).cursor - base_est(i).data_valid.cursor > 1e-3)      
end

%
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
end

%
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
end

%%
%Plot co-expression and BMI E1, E2 cells
%--------------------------------------------------------------------------
%Plot the red-green channel with bmi cells in blue
%Coexpression plot: 
% % save red and Im in folder/animal/day
% filetosave = fullfile(savePath, 'red.mat');
% save(filetosave,'Im', 'Img', 'red', 'redGreen', 'holoMask', 'holoMaskRedGreen')

mask_rg_data = load(im_file_cell{1}); 
%Im - red channel
%Img - green channel
green_scale = 7.7; 
red_scale   = 6.7; 

[fov_im, E2_mask, E1_mask, bmi_mask, E1_overlay, E2_overlay, BMI_overlay] = ...
    plot_coexpression(mask_rg_data, tinfo, red_scale, green_scale, plot_dir);

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
    hit_idxs    = base(i).valid_hit_times; 
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
win_zoom    = [5 10]; 
save_dir    = plot_dir; 
save_bool   = 0; 
plot_big    = 0; 
plot_block_neural(block_str, block_data, event_struct, frame_rate, win_zoom, save_dir, save_bool, plot_big)

%%
%BMI: 
%--------------------------------------------------------------------------
frame_rate  = 30; 
time_scale  = 1/(frame_rate*60); 
block_str   = 'bmi'; 
block_data  = bmi_est; 

%Define Events to Lock to: 

clear event_struct
for i = 1:length(block_data)
    valid_idxs                  = block_data(i).data_valid.valid_idxs;
    
    event_i = 0; 
    %Target Hit: 
    event_i     = event_i +1;
    
    hit_idxs    = find(bmi(i).data.selfHits(valid_idxs)); 
    hit_time                    = time_scale*hit_idxs; 
    
    event_struct(i,event_i).data      = hit_time;     
    event_struct(i,event_i).label     = 'hit'; %hit_b2b'; 
    event_struct(i,event_i).valid     = 1;     
    event_struct(i,event_i).line      = '--k';
end
bmi_event = event_struct; 

%%
win_zoom    = [0 5]; 
save_dir    = plot_dir; 
save_bool   = 0; 
plot_big    = 0; 
plot_block_neural(block_str, block_data, event_struct, frame_rate, win_zoom, save_dir, save_bool, plot_big)

%--------------------------------------------------------------------------
%%
%MANUALLY Make pretty plots of example data: 

save_bool = 0
win_zoom = [10 15]; 
%good examples: [10 15]

%Idxs to plot: 
% n_valid = find(~isnan(n_pad.pretrain(1,:)));
% start_offset_idx = n_valid(1);
start_offset_idx = 0; 
start_idx = start_offset_idx+win_zoom(1)*frame_rate*60; 
stop_idx = start_idx+(win_zoom(2)-win_zoom(1))*frame_rate*60;
idx_zoom = [start_idx stop_idx]; 

%
%Collect the neural data to plot
for i = 1:num_BMI
    n_cell{i} = pretrain_est.dff_z(i,start_idx:stop_idx);
    t_cell{i} = linspace(win_zoom(1), win_zoom(2), length(n_cell{i})); 
end

num_neurons = length(n_cell); 
offset = 0;
offset_vec = [offset]; 
h = figure; hold on;
for i=1:num_neurons
    data_i = n_cell{i}; 
    data_i = data_i(:); 
    
    t_i = t_cell{i}; 
    %First trace: 
    y_plot = data_i(:,1);
    first_trace_offset = min(y_plot); 
    y_plot = y_plot-first_trace_offset;
    y_amp = max(y_plot); 
    if(i>1)
        offset = offset + y_amp;
        offset_vec = [offset_vec offset]; 
    end
    y_plot = y_plot-offset;
    plot(t_i, y_plot, 'Color', E_color{E_id(i)}); 
end

evt_line = pretrain_event(2).line;
evt_data = pretrain_event(2).data;
evt_data_plot = evt_data(evt_data>= win_zoom(1) & evt_data <= win_zoom(2));
vline(evt_data_plot, evt_line);
set(gca,'TickDir','out');
%SAVE: 
if save_bool
    export_fig(h, fullfile(plot_dir, '5min_stim_example.png')); 
    export_fig(h, fullfile(plot_dir, '5min_stim_example.eps')); 
end

%%
%MANUALLY Make pretty plots of example data: 


save_bool = 1
win_zoom = [4.0 5.5]; 
%good examples: [11.5 13]

%Idxs to plot: 
% n_valid = find(~isnan(n_pad.pretrain(1,:)));
% start_offset_idx = n_valid(1);
start_offset_idx = 0; 
start_idx = start_offset_idx+win_zoom(1)*frame_rate*60; 
stop_idx = start_idx+(win_zoom(2)-win_zoom(1))*frame_rate*60;
idx_zoom = [start_idx stop_idx]; 

%
%Collect the neural data to plot

for i = 1:num_BMI
    n_cell{i} = pretrain_est.dff_z(i,start_idx:stop_idx);
    t_cell{i} = linspace(win_zoom(1), win_zoom(2), length(n_cell{i})); 
end

num_neurons = length(n_cell); 
n_plot = [1 5]; 
offset = 0;
offset_vec = [offset]; 
h = figure; hold on;
for i=1:num_neurons
    if sum(n_plot==i) > 0
        data_i = n_cell{i}; 
        data_i = data_i(:); 

        t_i = t_cell{i}; 
        %First trace: 
        y_plot = data_i(:,1);
        first_trace_offset = min(y_plot); 
        y_plot = y_plot-first_trace_offset;
        y_amp = max(y_plot); 
        if(i>1)
            offset = offset + y_amp;
            offset_vec = [offset_vec offset]; 
        end
        y_plot = y_plot-offset;
        plot(t_i, y_plot, 'Color', E_color{E_id(i)}); 
    end
end

evt_line = pretrain_event(2).line;
evt_data = pretrain_event(2).data;
evt_data_plot = evt_data(evt_data>= win_zoom(1) & evt_data <= win_zoom(2));
vline(evt_data_plot, evt_line);
set(gca,'TickDir','out');
%SAVE: 
if save_bool
    export_fig(h, fullfile(plot_dir, 'single_neuron_stim_example.png')); 
    export_fig(h, fullfile(plot_dir, 'single_neuron_stim_example.eps')); 
end

%%
%1) 
%base_event, pretrain_event, bmi_event
%

%%
%1) Let's see baseline, pretrain, bmi distributions of E2-E1
%2) Verification: reproduce calculated cursor with saved neural activity

%%
tinfo = load(target_info_cell{1}); 


%%
%1) Let's see baseline, pretrain, bmi distributions of E2-E1
base = load(cal_cell{1}); 
cursor.base = base.cursor_obs; 

%%
pretrain = load(pretrain_file_cell{1}); 
valid_idxs = find(~isnan(pretrain.data.cursor));
cursor.pretrain = pretrain.data.cursor(valid_idxs); 

%%
for i =1:length(bmi_file_cell)
    bmi(i) = load(bmi_file_cell{i}); 
end

%%
bmi_accum = []; 

bmi = load(bmi_file_cell{1}); 
c = bmi.data.cursor;
valid_idxs = find(~isnan(c)); 
c_v = c(valid_idxs); 
bmi_accum = [bmi_accum c_v];

bmi = load(bmi_file_cell{2}); 
c = bmi.data.cursor;
valid_idxs = find(~isnan(c)); 
c_v = c(valid_idxs); 
bmi_accum = [bmi_accum c_v];
cursor.bmi = bmi_accum; 

%%
min_data = min([cursor.base(:); cursor.pretrain(:); cursor.bmi(:)]); 
max_data = max([cursor.base(:); cursor.pretrain(:); cursor.bmi(:)]); 
num_bins = 100; 
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
h = figure;
hold on; 
plot(bin_centers, base_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
% plot(bin_centers, pretrain_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
% plot(bin_centers, bmi_bin_n, '.-', 'LineWidth', 1.5, 'MarkerSize', 10); 
vline(base.T, 'r:'); 
% legend({'base', 'pretrain', 'bmi', 'target'}); 

xlim([-1.5 1]); 
xlabel('cursor value: E2-E1'); 
ylabel('frac time entered'); 
%TODO: get hist code I used when doing fernandoBMI analysis
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, 'base_hist.eps')); 
%%
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
export_fig(h, fullfile(plot_dir, 'base_pretrain_bmi_hist.eps')); 
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
%Plot E1, E2 dff_smooth histograms

%baseline is easy, we have the dff smooth: 
n.base.E1 = base.n_analyze(:, find(base.E1_sel)); 
n.base.E2 = base.n_analyze(:, find(base.E2_sel)); 





%%
% %%
% %BUG: 
% %Pretrain:
% movingAverageFrames = 4; 
% cursor_zscore_bool = 0; 
% ind = 1; 
% pretrain = load(pretrain_file_cell{ind});
% [pretrain_est_bug] = ...
%     base_bmiAct2cursor(pretrain.data.baseVector, pretrain.data.bmiAct, movingAverageFrames, cursor_zscore_bool, pretrain.bData, pretrain.data.cursor);
% sum(pretrain_est.cursor - pretrain_est.data_valid.cursor > 1e-3)
% 

% %%
% %NOTBUG:
% %Pretrain:
% movingAverageFrames = 4; 
% cursor_zscore_bool = 0; 
% ind = 1; 
% pretrain = load(pretrain_file_cell{ind});
% [pretrain_est_nb] = ...
%     base_bmiAct2cursor(pretrain.data.baseVector, pretrain.data.bmiAct, movingAverageFrames, cursor_zscore_bool, pretrain.bData, pretrain.data.cursor);
% sum(pretrain_est.cursor - pretrain_est.data_valid.cursor > 1e-3)

%%
h = figure;
hold on; 
plot(pretrain_est_nb.dff(1,:)); 
plot(pretrain_est_bug.dff(1,:)); 
legend({'not bug', 'bug'}); 

%%
h = figure; 
hold on; 
plot(pretrain_est_nb.data_valid.baseVector(1,:)); 
plot(pretrain_est_nb.data_valid.baseVector(2,:)); 
plot(pretrain_est_nb.data_valid.baseVector(3,:)); 
%%
h = figure; 
hold on;
plot(pretrain_est.data_valid.baseVector(1,:)); 
plot(pretrain_est.data_valid.baseVector(2,:)); 

%%
h = figure;
scatter(pretrain_est_nb.cursor, pretrain_est_bug.cursor)

%%
h = figure;
scatter(pretrain_est_nb.dff(1,:), pretrain_est_bug.dff(1,:))



%%
%Collect BMI neural timeseries: 
n.base          = base.n_analyze.'; %which should be dff
n.pretrain      = pretrain_est.dff; 
n.bmi           = [bmi_est_cell{1}.dff bmi_est_cell{2}.dff];

n.base_z        = base.dff_z.'; 
n.pretrain_z    = pretrain_est.dff_z; 
n.bmi_z         = [bmi_est_cell{1}.dff_z bmi_est_cell{2}.dff_z];

%%
%Hist: 
holo_valid = find(pretrain.data.holoDelivery(pretrain_est.data_valid.valid_idxs)); 
holo_win = [];

num_frames_post_stim = 30; 

for i=1:length(holo_valid)
    holo_win = [holo_win holo_valid(i):(holo_valid(i)+num_frames_post_stim)]; 
end


holo_cursor = pretrain_est.cursor(holo_win); 


%%
 [data_binned, holo_cursor_bin] = bin_data(bin_edges, holo_cursor);
holo_cursor_bin_n = holo_cursor_bin/sum(holo_cursor_bin); 

%%
h = figure;
hold on; 
plot(bin_centers, pretrain_bin_n, '.-k', 'LineWidth', 1.5, 'MarkerSize', 10); 
plot(bin_centers, holo_cursor_bin_n, '.-r', 'LineWidth', 1.5, 'MarkerSize', 10); 
xlim([-1 1]); 
vline(base.T)
legend({'pretrain E2-E1', 'stim E2-E1'}); 
xlabel('E2-E1'); 
ylabel('frac time entered'); 
set(gca,'TickDir','out');
% export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_cursor.eps')); 

%%
n.pretrain_holo = n.pretrain(:, holo_win); 

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


% [data_binned, base_E1_bin] = bin_data(bin_edges_E1, n.base_E1);
[data_binned, pretrain_E1_bin] = bin_data(bin_edges_E1, n.pretrain_E1);
[data_binned, holo_E1_bin] = bin_data(bin_edges_E1, n.pretrain_holo_E1);

base_E1_bin_n = base_E1_bin/sum(base_E1_bin); 
pretrain_E1_bin_n = pretrain_E1_bin/sum(pretrain_E1_bin); 
holo_E1_bin_n = holo_E1_bin/sum(holo_E1_bin); 
%
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
% export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_E1.eps')); 



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

%
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
% export_fig(h, fullfile(plot_dir, 'pretrain_vs_holo_E2.eps')); 

%%
base_mean       = mean(n.base,2)
pretrain_mean   = mean(n.pretrain, 2)
bmi_mean        = mean(n.bmi, 2)

% base_mean =
% 
%     0.0393
%     0.0208
%     0.0018
%     0.0059
%     0.0021
%     0.0126
%     0.0026
%     0.0194
% 
% 
% pretrain_mean =
% 
%     0.0046
%    -0.0252
%    -0.0176
%    -0.0125
%    -0.0058
%    -0.0142
%    -0.0244
%    -0.0079
% 
% 
% bmi_mean =
% 
%     0.1051
%     0.1042
%     0.1088
%     0.0254
%     0.0548
%     0.0672
%     0.0251
%     0.0283
    
%%
%bmi1_mean = 
% 
%     0.0763
%     0.0903
%     0.0715
%    -0.0016
%     0.0322
%     0.0535
%     0.0127
%    -0.0097

%bmi2_mean = 
%
%     0.1084
%     0.1058
%     0.1131
%     0.0286
%     0.0574
%     0.0688
%     0.0266
%     0.0327
    
%%
%cov: row is observation
base_cov = cov(n.base.');
pretrain_cov = cov(n.pretrain.');
bmi_cov = cov(n.bmi.');

%%
pool_cov = [base_cov(:); pretrain_cov(:); bmi_cov(:)];
cmax = max(pool_cov) 
cmin = min(pool_cov)
%%
im = base_cov; 
h = figure;
imagesc(im); colorbar;
caxis([cmin cmax]); 
axis equal
title('base'); 

im = pretrain_cov; 
h = figure;
imagesc(im); colorbar;
caxis([cmin cmax]); 
axis equal
title('pretrain'); 

im = bmi_cov; 
h = figure;
imagesc(im); colorbar;
caxis([cmin cmax]); 
axis equal
title('bmi'); 
%%
base_z_mean       = mean(n.base_z,2)
pretrain_z_mean   = mean(n.pretrain_z, 2)
bmi_z_mean        = mean(n.bmi_z, 2)


% base_z_mean =
% 
%    1.0e-15 *
% 
%     0.3899
%     0.7092
%     0.1431
%     0.3448
%     0.1986
%    -0.2024
%    -0.1872
%    -0.4429
% 
% 
% pretrain_z_mean =
% 
%    -0.1085
%    -0.0922
%    -0.0782
%    -0.0712
%    -0.0772
%    -0.1285
%    -0.2848
%    -0.1125
% 
% 
% bmi_z_mean =
% 
%     0.2058
%     0.1670
%     0.4307
%     0.0760
%     0.5155
%     0.2614
%     0.2368
%     0.0363

%%
%Histogram the neural activity 

%%
%Plot some example data from each epoch: 
num_samples_plot = 20000; 
%Plot colors: 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
E_id = [ones(1,4) 2*ones(1,4)]'; 
% offset = 0; 

n_plot = n.base_z(:,1:num_samples_plot).'; 

[h, offset_vec] = plot_E_activity(n_plot, E_id, E_color);
title('base'); 

n_plot = n.pretrain_z(:,1:num_samples_plot).'; 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(n_plot, E_id, E_color);
title('pretrain'); 

n_plot = n.bmi_z(:,1:num_samples_plot).'; 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(n_plot, E_id, E_color);
title('bmi'); 

%%
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 

% n_plot = zscore(n.base.'); 
n_plot = n.base.'; 
t = (1:size(n_plot,1))/(26*60); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color, 0);
title('base'); 

set(gca,'TickDir','out');
% export_fig(h, fullfile(plot_dir, 'base.eps')); 


%
n_plot = n.pretrain.'; 
t = (1:size(n_plot,1))/(30*60); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color, 0);
title('pretrain'); 

%
n_plot = n.bmi.'; 
t = (1:size(n_plot,1))/(30*60); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color, 0);
title('bmi'); 

%%
pretrain_f0_valid = pretrain.data.baseVector(:, pretrain_est.data_valid.valid_idxs);
t = (1:size(n_plot,1))/(30*60); 
n_plot = pretrain_f0_valid.'; 
h = figure; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color, 0);

%%
pretrain_f_valid = pretrain.data.bmiAct(:, pretrain_est.data_valid.valid_idxs);
t = (1:size(n_plot,1))/(30*60); 
n_plot = pretrain_f_valid.'; 
h = figure; 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color, 0);

%%
h = figure;
hold on;
plot(pretrain_f_valid(2,:))
plot(pretrain_f0_valid(2,:), 'LineWidth', 5)
legend({'f', 'f0'}); 

%%
h = figure;
plot(n.pretrain(2,:)); 
title('dff'); 

%%
h = figure; 
plot(pretrain_est.dff(2,:)); 
title('dff est'); 

%%
%cov: row is observation
base_cov = cov(n.base_z.');
pretrain_cov = cov(n.pretrain_z.');
bmi_cov = cov(n.bmi_z.');

pool_cov = [base_cov(:); pretrain_cov(:); bmi_cov(:)];
cmax = max(pool_cov) 
cmin = min(pool_cov)
%%
im = base_cov; 
h = figure;
imagesc(im); colorbar;
caxis([cmin cmax]); 
axis equal
title('base'); 

im = pretrain_cov; 
h = figure;
imagesc(im); colorbar;
caxis([cmin cmax]); 
axis equal
title('pretrain'); 

im = bmi_cov; 
h = figure;
imagesc(im); colorbar;
caxis([cmin cmax]); 
axis equal
title('bmi'); 
%%
%Load the pretrain CSV, make a PSTH locked to stim: 
% %Figure out which channel is the laser trigger
% time_idx    = 1; 
% laser_idx   = 8; 
% h = figure; 
% plot(voltageRec((1:50000), 1), voltageRec(1:50000, 8))

%%
%

h = figure;
plot(pretrain.data.holoDelivery)

%%
n_plot = n.pretrain.'; 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(n_plot, E_id, E_color);
title('pretrain'); 

hold on; 
vline(event_idxs); 

%%
n_plot = n.pretrain.'; 
E_id = [ones(1,4) 2*ones(1,4)]'; 
[h, offset_vec] = plot_E_activity(n_plot, E_id, E_color);
title('pretrain'); 

hold on; 
vline(event_idxs); 

%%
h = figure;
hold on; 
plot(pretrain_est.data_valid.cursor); 
plot(pretrain_est.cursor); 
legend({'saved', 'recon'}); 

%%
h = figure;
scatter(pretrain.data.cursor(~isnan(pretrain.data.cursor)), pretrain_est.cursor); 
%%
n_plot = zscore(n.pretrain.'); 
n_plot = n_plot(1:10000, :); 
E_id = [ones(1,4) 2*ones(1,4)]'; 
t = (1:size(n_plot,1))/(30*60); 
[h, offset_vec] = plot_E_activity(t, n_plot, E_id, E_color, 0);
title('pretrain'); 

hold on; 
vline(event_idxs/(30*60)); 
set(gca,'TickDir','out');
xlabel('time (min)'); 

%%
data_mat = n.pretrain.'; 
event_idxs = find(pretrain.data.holoDelivery(pretrain_est.data_valid.valid_idxs));
win = [-500 500]; 
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
% export_fig(h, fullfile(plot_dir, ['pretrain_stim_psth_win' num2str(win_len) '.eps'])); 

%%
%Zscore, locked to stim
data_mat = zscore(n.pretrain.'); 
event_idxs = find(pretrain.data.holoDelivery(pretrain_est.data_valid.valid_idxs));
win = [-100 100]; 
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
ylabel('deltaF over F'); 
title('PSTH of Ensemble Activity Locked to 2p stimulation'); 
set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, ['pretrain_stim_psth_Z_win' num2str(win_len) '.eps'])); 


%%
% PSTH locked to self-hits in pretrain and BMI

%%
%PSTH locked to self-hits in pretrain: 
%%
% Pretrain selfhit
data_mat = n.pretrain.'; 
event_idxs = find(pretrain.data.selfHits(pretrain_est.data_valid.valid_idxs));
%Remove consecutive hits:
length(event_idxs)
event_idxs(find(diff(event_idxs)==1)+1) = []; 
length(event_idxs)


win = [-100 100]; 
win_len = win(2)-win(1); 
[psth_mean, psth_sem, psth_mat] = calc_psth(data_mat, event_idxs, win);

%
t_plot = (win(1):win(2))/30; 
h = figure; hold on;
offset = 0; 
for i=1:num_neurons
    y_plot = psth_mean(:,i); 
    y_plot = y_plot-min(y_plot);
    y_amp = max(y_plot+psth_sem(:,i)); 
    offset = offset + y_amp + 0.1; 
    y_sem = psth_sem(:,i)-min(y_plot); 
    
    plot(t_plot, y_plot-offset, 'Color', E_color{(E_id(i))}); 
    errbar(t_plot, y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
end
vline(0)
% vline((psth_win(2)-psth_win(1))/2+1); 
xlabel('time (sec)')
ylabel('deltaF over F'); 
title('PSTH of Ensemble Activity Locked to Self Hit'); 
set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, ['pretrain_selfHit_psth_win' num2str(win_len) '.eps'])); 


%%
%Pretrain selfhit, Zscore
data_mat = zscore(n.pretrain.'); 
event_idxs = find(pretrain.data.selfHits(pretrain_est.data_valid.valid_idxs));
length(event_idxs)
event_idxs(find(diff(event_idxs)==1)+1) = []; 

win = [-100 100]; 
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
ylabel('deltaF over F'); 
title('PSTH of Ensemble Activity Locked to Self Hit'); 
set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, ['pretrain_selfHit_psth_Z_win' num2str(win_len) '.eps'])); 



%%
%Process hits, so I only include hits with one sample of b2base in between 
hit_idxs    = find(pretrain.data.selfHits(pretrain_est.data_valid.valid_idxs));
b2base_idxs = find(pretrain_est.data_valid.cursor <= bData.T1);
valid_hit_idxs = [hit_idxs(1)];
for i = 2:length(hit_idxs)
    if sum((b2base_idxs > hit_idxs(i-1)) & (b2base_idxs <= hit_idxs(i))) > 0
        valid_hit_idxs = [valid_hit_idxs hit_idxs(i)]; 
    end
end

%%
%Sliding window of self-hit rate in pretrain and BMI: 

%%
close all

%%
%Pretrain: sliding window of self-hit: 
frameRate = 30; 
win = 3*60*frameRate; %
rate_filt = ones(win, 1)/(win/(frameRate*60)); 
valid_hit_idxs
hit_binary = pretrain.data.selfHits(pretrain_est.data_valid.valid_idxs); 
hit_rate = conv(hit_binary, rate_filt, 'valid'); 

hit_b2base_binary = zeros(length(hit_binary), 1); 
hit_b2base_binary(valid_hit_idxs) = 1; 
hit_b2base_rate = conv(hit_b2base_binary, rate_filt, 'valid'); 

%%
%Pretrain: sliding window of self-hit: 
frameRate = 30; 
win = 3*60*frameRate; %
rate_filt = ones(win, 1)/(win/(frameRate*60)); 
valid_hit_idxs
hit_binary = pretrain.data.selfHits(pretrain_est.data_valid.valid_idxs); 
hit_rate = conv(hit_binary, rate_filt, 'valid'); 

hit_b2base_binary = zeros(length(hit_binary), 1); 
hit_b2base_binary(valid_hit_idxs) = 1; 
hit_b2base_rate = conv(hit_b2base_binary, rate_filt, 'valid'); 

%%
t = (1:length(hit_rate))/(30*60); 
h = figure;
plot(t, hit_rate); 
xlabel('time (min)'); 
ylabel('hit rate (hits per min)'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, ['pretrain_selfHit_sliding_win' num2str(win) '.eps'])); 

%
t = (1:length(hit_rate))/(30*60); 
h = figure;
plot(t, hit_b2base_rate); 
xlabel('time (min)'); 
ylabel('hit rate (hits per min) with b2base'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, ['pretrain_selfHit_b2base_sliding_win' num2str(win) '.eps'])); 

%
h = figure;
% hold on; 
ax = plotyy(t, hit_rate, t, hit_b2base_rate); 
% plot(t, hit_b2base_rate); 
legend({'hits', 'hits with return to base'}); 
xlabel('time (min)'); 
% ylabel('hit rate (hits per min)'); 
ylabel(ax(1),'hit rate (hits per min)')
ylabel(ax(2),'hit rate with return to base (hits per min)')
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, ['pretrain_selfHit_both_sliding_win' num2str(win) '.eps'])); 

%%
%BMI: Sliding window of self-hit:
frameRate = 30; 
win = 3*60*frameRate; %
rate_filt = ones(win, 1)/(win/(frameRate*60)); 

%Pool over both BMI files: 
ind = 1; 
bmi = load(bmi_file_cell{ind}); 
hit1_binary = bmi.data.selfHits(bmi_est_cell{ind}.data_valid.valid_idxs);
ind = 2; 
bmi = load(bmi_file_cell{ind}); 
hit2_binary = bmi.data.selfHits(bmi_est_cell{ind}.data_valid.valid_idxs);
hit_binary = [hit1_binary hit2_binary]; 

hit_rate = conv(hit_binary, rate_filt, 'valid'); 

% h = figure;
% plot(hit_binary); 

t = (1:length(hit_rate))/(30*60); 
h = figure;
plot(t, hit_rate); 
xlabel('time (min)'); 
ylabel('hit rate (hits per min)'); 
set(gca,'TickDir','out');
ylim([0 25]); 
% export_fig(h, fullfile(plot_dir, ['BMI_selfHit_sliding_win' num2str(win) '.eps'])); 

%%
%Bar plot of self-hit rate comparing base, pretrain, BMI

%%
bmi1 = load(bmi_file_cell{1}); 
bmi2 = load(bmi_file_cell{2}); 
%%
% hit rate: 
%0.5
%388/40
%308/40

%%
h = figure;
bar([0.5 388/40 308/40]); 
ylabel('hits per min'); 
set(gca,'TickDir','out');
export_fig(h, fullfile(plot_dir, 'hits_per_min_bar.eps')); 

%%
%Sliding window of self-hit



%Cum rate

%Correct for back-to-baseline


%%
%PSTH of neural locked to hit in BMI mode: 
%Analyze first BMI file: 
z_bool = 1; 

bmi = load(bmi_file_cell{1}); 
bmi_est = bmi_est_cell{1}; 

valid_idxs = bmi_est.data_valid.valid_idxs; 

% BMI selfhit
% event_idxs = find(bmi.data.selfHits(valid_idxs))+length(bmi_est_cell{1}.data_valid.valid_idxs);
event_idxs = find(bmi.data.selfHits(valid_idxs));
win = [-100 100]; 
win_len = win(2)-win(1); 

if(z_bool)
    data_mat = zscore(n.bmi.'); 
    im_name = ['BMI_psthZ_win' num2str(win_len) '.eps'];
else
    data_mat = n.bmi.'; 
    im_name = ['BMI_psth_win' num2str(win_len) '.eps'];
end
[psth_mean, psth_sem, psth_mat] = calc_psth(data_mat, event_idxs, win);

%PLOTTING:
t_plot = (win(1):win(2))/30; 
h = figure; hold on;
offset = 0; 
for i=1:num_neurons
    y_plot = psth_mean(:,i); 
    y_plot = y_plot-min(y_plot);
    y_amp = max(y_plot+psth_sem(:,i)); 
    offset = offset + y_amp + 0.1; 
    y_sem = psth_sem(:,i)-min(y_plot); 
    
    plot(t_plot, y_plot-offset, 'Color', E_color{(E_id(i))}); 
    errbar(t_plot, y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
end
vline(0)
% vline((psth_win(2)-psth_win(1))/2+1); 
xlabel('time (sec)')
ylabel('deltaF over F'); 
title('PSTH of Ensemble Activity Locked to BMI Hit'); 
set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, im_name)); 



%%
%Choosing cells plot

%%
%Raw activity in baseline, pretrain, BMI