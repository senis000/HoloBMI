function load_path = define_and_load_bmi_paths_debug(home_dir)
i = 0;

i=i+1;
load_path(i).path = fullfile(home_dir, 'roi_process');
load_path(i).gen_bool = 0; 

i = i+1; 
load_path(i).path = fullfile(home_dir, 'baseline_target_calibration');
load_path(i).gen_bool = 0; 

i = i+1; 
load_path(i).path = fullfile(home_dir, 'matlab_exchange'); 
load_path(i).gen_bool = 0; 

i = i+1; 
load_path(i).path = fullfile(home_dir, 'baseline_target_calibration, plot_util')
load_path(i).gen_bool = 1; 

% i = i+1; 
% load_path(i).path = 'G:\VivekNuria\Code\dex';
% load_path(i).gen_bool = 1; 
% 
% i = i+1; 
% load_path(i).path = 'G:\VivekNuria\Code\Kenichi';
% load_path(i).gen_bool = 1; 

num_paths = length(load_path); 
for i_path = 1:num_paths
    if(load_path(i_path).gen_bool)
        addpath(genpath(load_path(i_path).path)); 
    else
        addpath(load_path(i_path).path); 
    end
end