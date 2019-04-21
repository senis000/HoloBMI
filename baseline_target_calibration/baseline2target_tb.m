
%% 
%ToDo:
%baseline file: will just have a matrix of neural activity
%%
addpath(genpath('/Users/vivekathalye/Dropbox/Code/holobmi_git/HoloBMI/baseline_target_calibration/plot_util')); 

%%
%Debug parameters: 
% (n_f_file, E_id, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)

save_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/results'
mkdir(save_dir); 

target_on_cov_bool = 0; 

frame_rate = 10; 
frames_per_reward_range = [1.5*60*frame_rate 1*60*frame_rate];
prefix_win = 40; %number frames to remove from the start of imaging
f0_win_bool = 1;
f0_win = 2*60*frame_rate; %frames
dff_win_bool = 1; 
dff_win = 2; %frames

data_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline'; 
n_f_file = fullfile(data_dir, 'baseline_IntegrationRois_00001_1.csv'); 
%Load baseline data:
A = csvread(n_f_file, 1, 0);
ts = A(:,1); 
frame = A(:,2); 
n = A(:, 3:end); 
num_n = size(n,2); 
[n_var, I] = sort(var(n, 0, 1), 'descend');
n = n(:,I); 
n(:,find(isnan(n_var))) = [];
n_var(find(isnan(n_var))) = [];

direct = [1 2 4 5 6 17 27 39]; %Hand chosen for debugging
n_f_raw = n(:,direct).'; %direct neuron fluorescence

E1_base = [6 17 27 39];
E2_base = [1 2 4 5];
%Save an n_f_file to be used for testing: 
base_test_path = fullfile(data_dir, 'base_data_test.mat');
f_base = n;
save(base_test_path, 'f_base'); 

%Transposed to be num_neurons x num_samples, as expected for the input to
%this file
% E_id = [2 2 2 2 1 1 1 1]'; 

%%
% % components file
Acomp_file = fullfile('/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/181213', 'redcomp.mat'); 
exist(Acomp_file)
comp_data = load(Acomp_file)

AComp = repmat(comp_data.AComp, [1, 2]); 

Acomp_save = fullfile('/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/181213', 'redcomp_test.mat'); 
save(Acomp_save, 'AComp'); 


%%
n_f_file = base_test_path; 
Acomp_file = fullfile('/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/181213', 'redcomp_test.mat'); 
baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)

%%


%%

% %%
% AComp = full(comp_data.AComp);
% %%
% im_idx =2;
% testim = reshape(AComp(:,im_idx), 256, 256);
% figure;
% imagesc(testim); 
% colormap(parula)
% colorbar
% axis square
% 
% %%
% h = figure;
% plot(comp_data.CComp(2,:))
% 
% %%
% comp_data.com
% 
% %%
% testim = reshape(sum(full(comp_data.AComp), 2), 256, 256);
% figure;
% imagesc(testim); 
% colormap(parula)
% colorbar
% axis square

%%
d_file = fullfile('/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/results', 'target_info.mat'); 
d = load(d_file)
%%
% %Debug parameters: 
% % (n_f_file, E_id, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
% 
% save_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/results'
% mkdir(save_dir); 
% 
% target_on_cov_bool = 0; 
% 
% frame_rate = 10; 
% frames_per_reward_range = [1.5*60*frame_rate 1*60*frame_rate];
% prefix_win = 40; %number frames to remove from the start of imaging
% f0_win_bool = 1;
% f0_win = 2*60*frame_rate; %frames
% dff_win_bool = 1; 
% dff_win = 2; %frames
% 
% data_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline'; 
% n_f_file = fullfile(data_dir, 'baseline_IntegrationRois_00001_1.csv'); 
% %Load baseline data:
% A = csvread(n_f_file, 1, 0);
% ts = A(:,1); 
% frame = A(:,2); 
% n = A(:, 3:end); 
% num_n = size(n,2); 
% [n_var, I] = sort(var(n, 0, 1), 'descend');
% n = n(:,I); 
% n(:,find(isnan(n_var))) = [];
% n_var(find(isnan(n_var))) = [];
% 
% direct = [1 2 4 5 6 17 27 39]; %Hand chosen for debugging
% n_f_raw = n(:,direct).'; %direct neuron fluorescence
% %Transposed to be num_neurons x num_samples, as expected for the input to
% %this file
% E_id = [2 2 2 2 1 1 1 1]'; 

%%
