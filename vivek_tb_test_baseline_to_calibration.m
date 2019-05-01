
%%
data_dir = 'F:\VivekNuria\a1\d1';

%%
%save(savePath + "BaselineOnline.mat", 'baseActivity')

base_file = fullfile(data_dir, 'BaselineOnline.mat')
exist(base_file)

%%
%TODO
%Need to streamline getting the files from 'BaselineAcqnvsPrairie.m' to 
%baseline2target.m

frame_rate = 30

n_f_file = base_file;
exist(n_f_file)
Acomp_file = fullfile(data_dir, 'redcomp.mat');
exist(Acomp_file)
E1_base = 1:4;
E2_base = 5:8;
frames_per_reward_range = [1 1.5]*60*frame_rate
target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 2*60*frame_rate
dff_win_bool = 1
dff_win = 2 %frames

save_dir = data_dir 
%%
baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, ...
    target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)