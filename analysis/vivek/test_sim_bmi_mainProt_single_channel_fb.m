%test BMI simulation: 
save_dir    = 'D:\Dropbox\Data\holobmi_nuria\NVI21_D07'
exist(save_dir)

data_dir    = 'D:\DATA\HoloBMI_2ndround\191111\NVI21\D07'
base_file   = fullfile(data_dir, 'target_calibration_ALL_20191111T174037.mat')
exist(base_file)
bmi_file    = fullfile(data_dir, 'BMI_online191111T182347.mat') 
exist(bmi_file)
n_f_file    = fullfile(data_dir, 'BaselineOnline191111T154732.mat') 
exist(n_f_file)

target_info_file = fullfile(data_dir, 'BMI_target_info_20191111T174037.mat'); 
exist(target_info_file); 

%%
%Sanity check baseline and bmi match: 
base_data = load(base_file)
bmi_data = load(bmi_file)
[task_settings] = define_BMI_task_settings();
[fb_settings]   = define_fb_audio_settings();


disp('E1 idxs: ')
base_data.E1_base
bmi_data.bData.E1_base

disp('E2 idxs: ')
base_data.E2_base
bmi_data.bData.E2_base

target_info = load(target_info_file);

%%
%Re-visit the analysis code I did:
%D:\Dropbox\Code\holobmi_git\HoloBMI\analysis\vivek\results_112019.m
%didn't say much
%%
%-n_f_file - contains matrix, neural fluorescence from baseline file, num_samples X num_neurons_baseline 

%%
%Re-simulate baseline calibration: 
recal_bool = 0
if recal_bool
    E1_base = base_data.E1_base;
    E1_base = base_data.E1_base;


    ndata = load(n_f_file);
    num_base_samples = sum(~isnan(ndata.baseActivity(1,:))); 
    baseline_frameRate = num_base_samples/(15*60);
    A_file = fullfile(data_dir, 'roi_data.mat'); 
    %roi_data_file; %fullfile(savePath, 'red.mat'); 
    exist(A_file)
    onacid_bool = 0

    sec_per_reward_range = [120 90]; 
    % sec_per_reward_range = [10 5]

    frames_per_reward_range = sec_per_reward_range*baseline_frameRate;
    disp('Time (s) per reward range: '); 
    disp(sec_per_reward_range); 
    disp('Frames per reward range: '); 
    disp(frames_per_reward_range)
    % sec_per_reward_range must be higher than 80seconds (to keep the
    % occurence of artificial vs natural higher than 80% 

    E2mE1_prctile = 98; 
    target_on_cov_bool = 0
    prefix_win = 40
    f0_win_bool = 1
    frameRate = 29.989;
    f0_win = 2*60*ceil(frameRate)
    dff_win_bool = 1
    dff_win = 4

    reward_per_frame_range = 1./frames_per_reward_range

    cursor_zscore_bool = 0;
    f0_init_slide = 0; 

    %Re-simulate baseline: 
    resim_base_dir = fullfile(save_dir, 'recal_base'); 
    [target_info_path, target_cal_ALL_path, fb_cal] = baseline2target_vE1strict_fb(n_f_file, A_file, onacid_bool,  ...
        E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
        prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, resim_base_dir, ...
        cursor_zscore_bool, f0_init_slide, E2mE1_prctile, fb_settings)
end

%%
% target_info = load(target_info_path)

%%
%resim base
resim_base_dir = fullfile(save_dir, 'resim_base'); 
ndata = load(n_f_file);
n_data = ndata.baseActivity; 
[result_save_path] = sim_bmi_vE1strict_fb(n_data, task_settings, target_info, resim_base_dir)

%%
sim_base_data = load(result_save_path)

%%
disp('params from base calibration'); 
base_data.T
base_data.E1_thresh
base_data.E2_subord_thresh

disp('params from sim base'); 
sim_base_data.T
sim_base_data.E1_thresh
sim_base_data.E2_subord_thresh

%%
disp('compare base cal results and sim base')
sim_base_data.num_valid_hits
base_data.num_valid_hits

sim_base_data.num_c1
base_data.num_c1

sim_base_data.num_c2
base_data.num_c2

sim_base_data.num_c3
base_data.num_c3

%%
%Simulate with the BMI data now: 
n_data = bmi_data.data.bmiAct; 
sim_dir = fullfile(save_dir, 'bmi_sim'); 
[result_save_path] = sim_bmi_vE1strict_fb(n_data, task_settings, target_info, sim_dir)

sim_bmi_data = load(result_save_path);

%%
sim_bmi_data.num_valid_hits
bmi_data.data.selfTargetCounter

% sim_base_data.num_c1
% base_data.num_c1
% 
% sim_base_data.num_c2
% base_data.num_c2
% 
% sim_base_data.num_c3
% base_data.num_c3