
%%
%Choose BMI Neurons by Plotting Baseline Activity
%First load: 
% 1) BaselineOnline.mat
% 2) redcomp.mat
% plotNeuronsBaseline(baseActivity, CComp, YrA, 20);
%%
%save(savePath + "BaselineOnline.mat", 'baseActivity')
base_file = fullfile(savePath, 'BaselineOnline.mat')
exist(base_file)

%%
rc = load(fullfile(savePath, 'redcomp.mat')); 
r   = load(fullfile(savePath, 'red.mat'));



%%
%TODO
%Need to streamline getting the files from 'BaselineAcqnvsPrairie.m' to 
%baseline2target.m


% Acomp_file = %fullfile(savePath, 'redcomp.mat');
%%
n_f_file = base_file;
exist(n_f_file)
ndata = load(n_f_file);
num_base_samples = sum(~isnan(ndata.baseActivity(1,:))); 
baseline_frameRate = num_base_samples/(15*60);
A_file = fullfile(savePath, 'red.mat'); 
exist(A_file)
onacid_bool = 0

% todo plot the traces from the baseline and from the stim
E1_base = [15 16 17 25];
E2_base = [11 9 24 6]; %Needs to have chrome!
sec_per_reward_range = [100 80];  

% baseline_frameRate = 
frames_per_reward_range = sec_per_reward_range*baseline_frameRate %[1 1.5]*60*frameRate
%multiply by frames per minute to convert
%to 

target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 2*60*ceil(frameRate)
dff_win_bool = 1
dff_win = 2
 
reward_per_frame_range = 1./frames_per_reward_range
%%
close all
baseline2target(n_f_file, A_file, onacid_bool, E1_base, E2_base, frames_per_reward_range, ...
    target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath)