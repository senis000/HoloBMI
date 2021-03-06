%%
home_dir = '/Users/vivekathalye/Dropbox/Code/holobmi_git/HoloBMI'
cd(home_dir); 
load_path = define_and_load_bmi_paths_debug(home_dir)
matlab_exchange_home = '/Users/vivekathalye/Dropbox/Code/MatlabExchange'; 
addpath(genpath(matlab_exchange_home)); 
%%
animal          = 'NVI16'
home_dir = '/Users/vivekathalye/Dropbox/Data/holo_sep2019';
data_dir = fullfile(home_dir, animal); 
debug_path = fullfile(data_dir, 'debug_data')
exist(debug_path)

day             = 'D0'

roi_data_file   = fullfile(data_dir, 'roi_data.mat'); 
base_file       = fullfile(data_dir, 'BaselineOnline190925T204015.mat'); 
bmi_file        = fullfile(data_dir, 'BMI_online190925T214408.mat'); 
cal_file        = fullfile(data_dir, 'BMI_target_info_20190925T205822.mat');
cal_ALL_file    = fullfile(data_dir, 'target_calibration_ALL_20190925T205822.mat'); 

exist(base_file)
exist(bmi_file)
exist(cal_file)
exist(cal_ALL_file)

strcMask_file   = fullfile(data_dir, 'strcMask.mat');
mask_data       = load(strcMask_file)

E1_base         = mask_data.E_base_sel(mask_data.E_id==1); 
E2_base         = mask_data.E_base_sel(mask_data.E_id==2);
%
%Constants:
frameRate = 29.989;

%% Calibrate Target with Baseline simulation
%--------------------------------------------------------------------------
%D0: (nothing)
%1) Parameters: 
% - sec_per_reward_range
% - f0_win (F0: how many frames to average over)
% - dff_win (F for Dff: how many frames to average over)
%--------------------------------------------------------------------------

% base_file = fullfile(savePath, 'BaselineOnline190514T221822.mat')
% base_file = bmi_base;
savePath            = debug_path; 

exist(base_file)
n_f_file            = base_file;
ndata               = load(n_f_file);
num_base_samples    = sum(~isnan(ndata.baseActivity(1,:))); 
baseline_frameRate  = num_base_samples/(15*60);
A_file              = roi_data_file; %fullfile(savePath, 'red.mat'); 
exist(A_file)
onacid_bool         = 0

sec_per_reward_range = [120 100]; 


frames_per_reward_range = sec_per_reward_range*baseline_frameRate;
disp('Time (s) per reward range: '); 
disp(sec_per_reward_range); 
disp('Frames per reward range: '); 
disp(frames_per_reward_range)
% sec_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 

target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 2*60*ceil(frameRate)
dff_win_bool = 1
dff_win = 4
 
reward_per_frame_range = 1./frames_per_reward_range
cursor_zscore_bool = 0;
f0_init_slide = 0; 

%%
E2mE1_prctile = 98; 

rerun_cal_bool = 1; 
if rerun_cal_bool
    close all
    [target_info_path, target_cal_ALL_path] = baseline2target_vE1strict(n_f_file, A_file, onacid_bool,  ...
        E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
        prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
        cursor_zscore_bool, f0_init_slide, E2mE1_prctile);
end

%%
base_data = load(base_file); 
cal_data = load(cal_file)
cal_ALL_data = load(cal_ALL_file);
%%
test_base_bool = 1;
test_bmi_bool = 0; 

if(test_bmi_bool) 
    bmi_data        = load(bmi_file); 
    debug_input     = bmi_data.data.bmiAct; 
    valid_idxs      = find(~isnan(debug_input(1,:))); 
    last_valid_idx  = max(valid_idxs)
    debug_input     = debug_input(:,1:last_valid_idx); 
    size(debug_input)
end

if(test_base_bool)
    debug_input = ...
        base_data.baseActivity(mask_data.E_base_sel,:);
    debug_input = debug_input(:, ~isnan(debug_input(1,:)));
    size(debug_input)
end


%%
% bmi_data.data.selfTargetCounter

%%
cal_ALL_data.num_valid_hits

%%

plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
t_plot = 1:size(debug_input,2);
[h, offset_vec] = plot_E_activity(t_plot.', debug_input.', mask_data.E_id, E_color, 0);

% [h, offset_vec] = plot_E_activity(t,n, E_id, E_color, offset)


%%
sim_bmi_bool = 1
single_bool         = 1; 
if single_bool
    if(test_bmi_bool)
        folder = fullfile(home_dir, animal, 'bmi_sim_single')
        mkdir(folder)
    else
%         folder = '/Users/vivekathalye/Dropbox/Data/holo_sep2019/NVI12/debug_data'
        folder = fullfile(home_dir, animal, 'debug_data'); 
    end
    exist(folder)
else
    if(test_bmi_bool)
        folder = fullfile(home_dir, animal, 'bmi_sim_double')
        mkdir(folder)
    else
        folder = fullfile(home_dir, animal, 'debug_data_v2'); 
%         '/Users/vivekathalye/Dropbox/Data/holo_sep2019/NVI12/debug_data_v2'
    end
    exist(folder)
end

% %%
% drop_prefix_bool = 0
% if drop_prefix_bool
%     bmi_sim_input = debug_input(:,(prefix_win+1):end); 
% else
%     bmi_sim_input = debug_input; 
% end
% 
% if sim_bmi_bool
%     expt_str = 'BMI'; 
%     baselineCalibrationFile = cal_file;
%     vectorHolo = []; 
%     vectorVTA = []; 
%     cursor_zscore_bool = 0;
%     debug_bool = 1; 
%     baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
% 
% %     BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
% %         expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
% %         cursor_zscore_bool, debug_bool, bmi_sim_input, baseValSeed, single_bool);    
%      BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4_OLD_TEST(folder, animal, day, ...
%         expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%         cursor_zscore_bool, debug_bool, debug_input, baseValSeed)    
% end


%%
single_bool         = 1; 
if single_bool
    if(test_bmi_bool)
        folder = fullfile(home_dir, animal, 'bmi_sim_single')
        mkdir(folder)
    end
    exist(folder)
else
    if(test_bmi_bool)
        folder = fullfile(home_dir, animal, 'bmi_sim_double')
        mkdir(folder)
    end
    exist(folder)
end

drop_prefix_bool = 0
if drop_prefix_bool
    bmi_sim_input = debug_input(:,(prefix_win+1):end); 
else
    bmi_sim_input = debug_input; 
end


expt_str = 'BMI'; 
baselineCalibrationFile = cal_file;
vectorHolo = []; 
vectorVTA = []; 
cursor_zscore_bool = 0;
debug_bool = 1; 
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan

BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, bmi_sim_input, baseValSeed);

%%
bmi_dir = fullfile(folder, animal, day)
% bmi_dir = '/Users/vivekathalye/Dropbox/Data/holo_sep2019/NVI16/bmi_sim_single/NVI16/D0'
% bmi_path = fullfile(bmi_dir, 'BMI_online190927T101911.mat'); 
% bmi_path = fullfile(bmi_dir, 'BMI_online190926T205945.mat'); 


%BMIAcqn first test: 
bmi_path = fullfile(bmi_dir, 'BMI_online190926T205945.mat'); 
exist(bmi_path)
bmi_sim0 = load(bmi_path)
bmi_sim0.data.selfTargetCounter


%BMIAcqn OLD TEST, to confirm first test:
% bmi_path1 = fullfile(bmi_dir, 'BMI_online190928T005845.mat'); 
% bmi_sim1 = load(bmi_path)
% bmi_sim1.data.selfTargetCounter
% NOTE: this one worked fine, which confused me.  It could be the
% dff2cursor edit that improved the code.  
%

bmi_path1 = fullfile(bmi_dir, 'BMI_online190928T095824.mat'); 
bmi_sim1 = load(bmi_path)
bmi_sim1.data.selfTargetCounter

%%

h = figure;
hold on;
% plot(cal_ALL_data.cursor_obs, 'LineWidth', 2); 
plot(bmi_sim0.data.cursor, 'LineWidth', 2); 
plot(bmi_sim1.data.cursor); 


%%
h = figure;
hold on;
plot(bmi_data.data.cursor, 'LineWidth', 2); 
plot(bmi_sim0.data.cursor); 
legend('bmi', 'sim'); 

%%
h = figure;
hold on;
plot(bmi_data.data.cursor);  
hline(bmi_data.bData.T1); 
plot(bmi_sim0.data.selfHits, 'k', 'LineWidth', 2); 
plot(bmi_data.data.selfVTA);

%%
h = figure;
hold on; 
plot(bmi_sim0.data.cursor); 
plot(1.5*bmi_sim0.data.c1_bool); 
plot(1.5*bmi_sim0.data.c2_bool); 
plot(1.5*bmi_sim0.data.c3_bool); 
plot(bmi_sim0.data.selfVTA, 'k', 'LineWidth', 2); 
legend('cursor' , 'c1', 'c2', 'c3', 'target'); 

%%
h = figure;
hold on; 
plot(bmi_sim0.data.c2_val); 
plot(1.5*bmi_sim0.data.c2_bool); 
plot(1.5*bmi_sim0.data.c1_bool); 
plot(bmi_sim0.data.selfVTA, 'k', 'LineWidth', 2); 
hline(bmi_sim0.bData.E1_thresh); 

%%
valid_hit_idxs = find(bmi_sim0.data.selfVTA)

%%
ind = 10
bmi_sim0.data.c1_bool(valid_hit_idxs(ind))
bmi_sim0.data.c2_bool(valid_hit_idxs(ind))
bmi_sim0.data.c3_bool(valid_hit_idxs(ind))
%%
bmi_dir = '/Users/vivekathalye/Dropbox/Data/holo_sep2019/NVI16/debug_data_v2/NVI16/D0'
bmi_path = fullfile(bmi_dir, 'BMI_online190927T115234.mat'); 
exist(bmi_path)
bmi_sim1 = load(bmi_path)
bmi_sim1.data.selfTargetCounter

%%
%sanity checks on SIM BMI with BASELINE DATA: 
%check 'data.bmiAct' is equal to debut_input: 
size(debug_input)
size(bmi_sim.data.bmiAct(:,1:length(debug_input)))
% find(debug_input ~= bmi_sim.data.bmiAct(:,1:length(debug_input)))
%%
cal_num_valid = length(debug_input)
bmi_num_valid = sum(~isnan(bmi_sim.data.bmiAct(1,:)))

%%
bmi_valid_idxs = find(~isnan(bmi_sim.data.bmiAct(1,:))); 
bmi_last_valid = bmi_valid_idxs(end);
%%
%NOTE: debug_input: num_neurons x num_samples
n_sel = 1; 
h = figure
hold on;
plot(debug_input(n_sel, :)); 
plot(bmi_sim.data.bmiAct(n_sel, 1:bmi_last_valid)); 
legend('base', 'bmi'); 
%in bmi_sim:
% NaN in the beginning for prefix_win number of samples

%%
%length(f0) = length(f_raw) - prefix_win*30
%Because of 'prefix_win', 
% Why is f0 length so much less than f_raw ???

% size(cal_ALL_data.f0)
% 
% ans =
% 
%        21510           8
% 
% size(cal_ALL_data.f_raw)
% 
% ans =
% 
%        25109           8

%%
%NOTE: 
%bmi_sim.data.baseVector: num_neurons X num_samples
%cal_ALL_data.f0: num_samples X num_neurons�
%
% t_cal = 1:cal_num_valid;
% t_bmi = 40+(1:bmi_num_valid); 

n_sel = 4; 

t_bmi = prefix_win +(f0_win:bmi_num_valid);
t_cal = 1:length(cal_ALL_data.f0);

cal_plot = cal_ALL_data.f0(t_cal,n_sel);
bmi_plot = bmi_sim.data.baseVector(n_sel,t_bmi);

h = figure;
hold on;
plot(cal_plot); 
plot(bmi_plot); 
legend({'cal', 'sim'})

%%
%F Gross comparison:
y_plot = cal_ALL_data.f0; 
t_plot = 1:size(y_plot,1);
[h, offset_vec] = plot_E_activity(t_plot, y_plot, mask_data.E_id, E_color, 0);
title('base cal'); 

t_bmi = prefix_win +(f0_win:bmi_num_valid);
y_plot = bmi_sim.data.baseVector(:,t_bmi).';
[h, offset_vec] = plot_E_activity(t_bmi, y_plot, mask_data.E_id, E_color, 0);
title('bmi sim'); 
% function [h, offset_vec] = plot_E_activity(t,n, E_id, E_color, offset)
% %4.18.19
% %n - neural activity, num_samples x num_neurons

%%
%Fsmooth comparison:

%NOTE: fsmooth is being calculated from sample 41 to 3639. at 3640, cursor becomes
%valid
h = figure;
hold on;
plot(bmi_sim.data.fsmooth(1,:)); 
plot(bmi_sim.data.cursor); 

%%
%Fsmooth comparison:

size(cal_ALL_data.f_smooth)
size(bmi_sim.data.fsmooth)
%%
%Fsmooth comparison:
n_sel = 1; 

% bmi_num_valid = sum(~isnan(bmi_sim.data.fsmooth(1,:))) 
% t_bmi = prefix_win +(f0_win:bmi_num_valid);
% t_cal = 1:length(cal_ALL_data.f0);

bmi_offset = 3600+2; 
cal_plot = cal_ALL_data.f_smooth(:,n_sel);
bmi_plot = bmi_sim.data.fsmooth(n_sel,bmi_offset:end);

h = figure;
hold on;
plot(cal_plot, 'LineWidth',2); 
plot(bmi_plot); 
legend({'cal', 'sim'})
title(['fsmooth: ' num2str(n_sel)]); 

%%
%dff comparison:
n_sel = 1; 

% bmi_num_valid = sum(~isnan(bmi_sim.data.fsmooth(1,:))) 
% t_bmi = prefix_win +(f0_win:bmi_num_valid);
% t_cal = 1:length(cal_ALL_data.f0);

bmi_offset = 3600+2; 
cal_plot = cal_ALL_data.dff(:,n_sel);
bmi_plot = bmi_sim.data.dff(n_sel,bmi_offset:end);

h = figure;
hold on;
plot(cal_plot, 'LineWidth',2); 
plot(bmi_plot); 
legend({'cal', 'sim'})
title(['DFF: ' num2str(n_sel)]); 

%%
%E1 comparison:

% bmi_num_valid = sum(~isnan(bmi_sim.data.fsmooth(1,:))) 
% t_bmi = prefix_win +(f0_win:bmi_num_valid);
% t_cal = 1:length(cal_ALL_data.f0);

c2_cal = zeros(size(cal_ALL_data.cursor_obs)); 
c2_cal(cal_ALL_data.c2)=1; 


bmi_offset = 3600+2; 
cal_plot = cal_ALL_data.E1_mean_analyze;
bmi_plot = bmi_sim.data.c2_val(n_sel,bmi_offset:end);

h = figure;
hold on;
plot(cal_plot, 'k'); 
plot(bmi_plot, 'r'); 
plot(c2_cal, 'k', 'LineWidth',2); 
plot(bmi_sim.data.c2_bool(n_sel,bmi_offset:end), 'r'); 
% legend({'cal', 'sim'})
title(['E1']); 

%%
%E1 bool:
h = figure;
hold on;
plot(c2_cal, 'k', 'LineWidth',2); 
plot(bmi_sim.data.c2_bool(1,bmi_offset:end), 'r'); 


%%
%E2 comparison:

% bmi_num_valid = sum(~isnan(bmi_sim.data.fsmooth(1,:))) 
% t_bmi = prefix_win +(f0_win:bmi_num_valid);
% t_cal = 1:length(cal_ALL_data.f0);

c3_cal = zeros(size(cal_ALL_data.cursor_obs)); 
c3_cal(cal_ALL_data.c3)=1; 


bmi_offset = 3600+2; 
cal_plot = cal_ALL_data.E2_subord_mean_analyze;
bmi_plot = bmi_sim.data.c3_val(n_sel,bmi_offset:end);

h = figure;
hold on;
plot(cal_plot, 'k'); 
plot(bmi_plot, 'r'); 
% plot(c2_cal, 'k', 'LineWidth',2); 
% plot(bmi_sim.data.c2_bool(n_sel,bmi_offset:end), 'r'); 
% legend({'cal', 'sim'})
title(['E2']); 


%%
%E2 bool: 
c3_cal = zeros(size(cal_ALL_data.cursor_obs)); 
c3_cal(cal_ALL_data.c3)=1; 

h = figure;
hold on;
plot(c3_cal, 'k', 'LineWidth',2); 
plot(bmi_sim.data.c3_bool(1,bmi_offset:end), 'r', 'LineWidth', 1.2); 


%%
%Cursor comparison: 
%cal_ALL_data
%bmi_sim.data
disp('base cal: ')
size(cal_ALL_data.cursor_obs)
% 
% disp('bmi sim: ')
% size(cal_ALL_data.cursor_obs)


h = figure;
hold on;
plot(cal_ALL_data.cursor_obs); 

% t_bmi = prefix_win +(f0_win:bmi_num_valid);
bmi_valid_idxs = find(~isnan(bmi_sim.data.cursor)); 
plot(bmi_sim.data.cursor(bmi_valid_idxs)); 
legend({'base cal', 'bmi sim'}); 

%%
%Confirm calibration calculates target_hits in equivalent way to BMI: 
size(cal_ALL_data.n_analyze)

% [dff_z, cursor, target_hit, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
%     dff2cursor_target(dff, bData, cursor_zscore_bool) 

bData               = cal_data;
cursor_zscore_bool  = 0; 
test_win            = [2000 4000];
win_samples         = test_win(1):test_win(2); 
num_samples_test    = test_win(2)-test_win(1)+1; 

%Reminder: 
% c1: cursor
% c2: E1_mean > E1_thresh
% c3: E2_subord_mean > E2_subord_thresh

test_data.cursor        = ones(1, num_samples_test);
test_data.target_hit    = ones(1, num_samples_test); 
test_data.c1_bool       = ones(1, num_samples_test);  
test_data.c2_bool       = ones(1, num_samples_test); 
test_data.c2_val        = ones(1, num_samples_test); 
test_data.c3_bool       = ones(1, num_samples_test); 
test_data.c3_val        = ones(1, num_samples_test); 
for i = 1:num_samples_test
    s_i                     = win_samples(i); 
    dff_i                   = cal_ALL_data.n_analyze(s_i,:); 
	[~, cursor, target_hit, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
        dff2cursor_target(dff_i, bData, cursor_zscore_bool);
    %ASSIGN:
    test_data.cursor(i)        = cursor;
    test_data.target_hit(i)    = target_hit;
    test_data.c1_bool(i)       = c1_bool;
    test_data.c2_bool(i)       = c2_bool;
    test_data.c2_val(i)        = c2_val;
    test_data.c3_bool(i)       = c3_bool;
    test_data.c3_val(i)        = c3_val;  
end

%%
%Cursor confirmation:

h = figure;
hold on; 
plot(cal_ALL_data.cursor_obs(win_samples), 'LineWidth', 2); 
plot(test_data.cursor); 
legend('base cal', 'confirm');

%%
%E1 confirmation:

c2_cal = zeros(size(cal_ALL_data.cursor_obs)); 
c2_cal(cal_ALL_data.c2)=1; 
c2_cal_win = c2_cal(win_samples); 

cal_plot    = cal_ALL_data.E1_mean_analyze(win_samples);
test_plot   = test_data.c2_val;

h = figure;
hold on; 
plot(c2_cal_win, 'k', 'LineWidth',2); 
plot(test_data.c2_bool, 'r'); 
plot(cal_plot, 'k', 'LineWidth',2); 
plot(test_plot, 'r'); 

%%
%E2 confirmation:

c3_cal = zeros(size(cal_ALL_data.cursor_obs)); 
c3_cal(cal_ALL_data.c3)=1; 
c3_cal_win = c3_cal(win_samples); 

cal_plot    = cal_ALL_data.E2_subord_mean_analyze(win_samples);
test_plot   = test_data.c3_val;

h = figure;
hold on; 
plot(c3_cal_win, 'k', 'LineWidth',2); 
plot(test_data.c3_bool, 'r'); 
plot(cal_plot, 'k', 'LineWidth',2); 
plot(test_plot, 'r'); 

%%
%CAL: Plot target hits 
h = figure;
hold on; 
vline(cal_ALL_data.hit_times, 'k')
vline(cal_ALL_data.valid_hit_times)
plot(cal_ALL_data.cursor_obs, 'LineWidth', 2);

%%
%BMI SIM: Plot target hits 
bmi_valid_idxs = find(~isnan(bmi_sim.data.cursor)); 
bmi_hits = find(bmi_sim.data.selfHits(bmi_valid_idxs)); 

h = figure;
hold on; 
% vline(bmi_sim.data.hit_times, 'k')
vline(bmi_hits);
plot(bmi_sim.data.cursor(bmi_valid_idxs), 'LineWidth', 2);

%%
%CAL + BMI: 
h = figure;
hold on; 
vline(cal_ALL_data.valid_hit_times, 'k')
vline(bmi_hits);
plot(cal_ALL_data.cursor_obs, 'LineWidth', 2);
plot(bmi_sim.data.cursor(bmi_valid_idxs), 'LineWidth', 2);

%%
%ToDo: add data to be saved if debugging into bmi code






































