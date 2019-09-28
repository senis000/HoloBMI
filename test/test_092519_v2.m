%%
home_dir = '/Users/vivekathalye/Dropbox/Code/holobmi_git/HoloBMI'
cd(home_dir); 
load_path = define_and_load_bmi_paths_debug(home_dir)
matlab_exchange_home = '/Users/vivekathalye/Dropbox/Code/MatlabExchange'; 
addpath(genpath(matlab_exchange_home)); 
%%
data_dir = '/Users/vivekathalye/Dropbox/Data/holo_sep2019'
debug_path = fullfile(data_dir, 'debug_data_v2')
exist(debug_path)

animal          = 'NVI16'
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

close all
[target_info_path, target_cal_ALL_path] = baseline2target_vE1strict(n_f_file, A_file, onacid_bool,  ...
    E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
    prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
    cursor_zscore_bool, f0_init_slide);

%%
base_data = load(base_file); 
cal_data = load(cal_file)
cal_ALL_data = load(cal_ALL_file); 

debug_input = ...
    base_data.baseActivity(mask_data.E_base_sel,:);
debug_input = debug_input(:, ~isnan(debug_input(1,:)));
size(debug_input)

%%

plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
t_plot = 1:size(debug_input,2);
[h, offset_vec] = plot_E_activity(t_plot.', debug_input.', mask_data.E_id, E_color, 0);

% [h, offset_vec] = plot_E_activity(t,n, E_id, E_color, offset)

%%
% folder = debug_path; 
folder = fullfile('/Users/vivekathalye/Dropbox/Data/holo_sep2019', 'debug_data_v2');
expt_str = 'BMI'; 
baselineCalibrationFile = cal_file;
vectorHolo = []; 
vectorVTA = []; 
cursor_zscore_bool = 0;
debug_bool = 1; 
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan

bmi_sim_input = debug_input(:,(prefix_win+1):end); 
BMIAcqnvsPrairienoTrialsHoloCL_debug_enable_v4(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, bmi_sim_input, baseValSeed);

%%
bmi_dir = '/Users/vivekathalye/Dropbox/Data/holo_sep2019/debug_data_v2/NVI16/D0'
bmi_path = fullfile(bmi_dir, 'BMI_online190927T115234.mat'); 
exist(bmi_path)
%%
bmi_sim = load(bmi_path)

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
%NOTE: debug_input: num_neurons x num_samples
n_sel = 1; 
h = figure
hold on;
plot(debug_input(n_sel, :)); 
plot(bmi_sim.data.bmiAct(n_sel, 1:bmi_num_valid)); 
legend('base', 'bmi'); 

%%
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
       
       


%Difference is due to 'prefix_win'

%%
%NOTE: 
%bmi_sim.data.baseVector: num_neurons X num_samples
%cal_ALL_data.f0: num_samples X num_neuronsﬂ
%
% t_cal = 1:cal_num_valid;
% t_bmi = 40+(1:bmi_num_valid); 

n_sel = 3; 

t_bmi = (prefix_win+f0_win):bmi_num_valid;
t_cal = 1:length(cal_ALL_data.f0);

cal_plot = cal_ALL_data.f0(t_cal,n_sel);
bmi_plot = bmi_sim.data.baseVector(n_sel,t_bmi);


h = figure;
hold on;
plot(cal_plot); 
plot(bmi_plot); 
legend({'cal', 'sim'})

%%
%Check the baseline fluorescence is similar to calibration

%%
bmi_true = load(bmi_file); 






