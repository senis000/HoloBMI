
%%
%Choose BMI Neurons by Plotting Baseline Activity
%First load: 
% 1) BaselineOnline.mat
% 2) redcomp.mat
plotNeuronsBaseline(baseActivity, CComp, YrA, 20);
%%
%save(savePath + "BaselineOnline.mat", 'baseActivity')
base_file = fullfile(savePath, 'BaselineOnline.mat')
exist(base_file)

%%
%TODO
%Need to streamline getting the files from 'BaselineAcqnvsPrairie.m' to 
%baseline2target.m

n_f_file = base_file;
exist(n_f_file)
Acomp_file = fullfile(savePath, 'redcomp.mat');
exist(Acomp_file)
E1_base = [27 23];
E2_base = [21 8]; %Needs to have chrome!
frames_per_reward_range = [1.5 0.1]*60*frameRate %[1 1.5]*60*frameRate

target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 10 %2*60*frameRate
dff_win_bool = 1
dff_win = 2 %frames
 
reward_per_frame_range = 1./frames_per_reward_range
%%
close all
baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, ...
    target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath)