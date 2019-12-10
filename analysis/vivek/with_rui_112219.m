
%%
%Cortex: 
clear n_base n_bmi
data_dir = '/Users/vivekathalye/Dropbox/Data/bmi_test_nov2019/fb_performance_check_112219/NVI18' 
% cal_file        = fullfile(data_dir, 'BMI_cal_20191119T095951.mat');
cal_ALL_file    = fullfile(data_dir, 'BMI_cal_ALL_20191122T100403.mat');
bmi_file        = fullfile(data_dir, 'BMI_online191122T104206.mat');
base_file       = fullfile(data_dir, 'BaselineOnline191122T094418.mat'); 
strcMask_file   = fullfile(data_dir, 'BMI_roi.mat');

exist(cal_ALL_file)
exist(bmi_file)
exist(strcMask_file)

cal_ALL_data    = load(cal_ALL_file); 
bmi_data        = load(bmi_file); 
base_data       = load(base_file); 
mask_data       = load(strcMask_file)


%%
%112219
clear n_base n_bmi
data_dir = '/Users/vivekathalye/Dropbox/Data/bmi_test_nov2019/fb_performance_check_112219/NY127/2019-11-22' 
% cal_file        = fullfile(data_dir, 'BMI_cal_20191119T095951.mat');
cal_ALL_file    = fullfile(data_dir, 'BMI_cal_ALL_20191122T112510.mat');
bmi_file        = fullfile(data_dir, 'BMI_online191122T120838.mat');
base_file       = fullfile(data_dir, 'BaselineOnline191122T110704.mat'); 
strcMask_file   = fullfile(data_dir, 'BMI_roi.mat');

exist(cal_ALL_file)
exist(bmi_file)
exist(strcMask_file)

cal_ALL_data    = load(cal_ALL_file); 
bmi_data        = load(bmi_file); 
base_data       = load(base_file); 
mask_data       = load(strcMask_file)


%%
%112119
%This works: 
clear n_base n_bmi
data_dir = '/Users/vivekathalye/Dropbox/Data/bmi_test_nov2019/fb_performance_check_112219/NY127/2019-11-21' 
% cal_file        = fullfile(data_dir, 'BMI_cal_20191119T095951.mat');

base_file       = fullfile(data_dir, 'BaselineOnline191121T094126.mat'); 
cal_ALL_file    = fullfile(data_dir, 'BMI_cal_ALL_20191121T095838.mat');
bmi_file        = fullfile(data_dir, 'BMI_online191121T102019.mat');

strcMask_file   = fullfile(data_dir, 'BMI_roi.mat');

exist(cal_ALL_file)
exist(bmi_file)
exist(strcMask_file)

cal_ALL_data = load(cal_ALL_file); 
bmi_data = load(bmi_file); 
base_data = load(base_file); 
mask_data       = load(strcMask_file)

%%
%111919
%This works: 
clear n_base n_bmi
data_dir = '/Users/vivekathalye/Dropbox/Data/bmi_test_nov2019/fb_performance_check_112219/NY127/2019-11-19' 
% cal_file        = fullfile(data_dir, 'BMI_cal_20191119T095951.mat');

base_file       = fullfile(data_dir, 'BaselineOnline191119T094038.mat'); 
cal_ALL_file    = fullfile(data_dir, 'BMI_cal_ALL_20191119T095951.mat');
bmi_file        = fullfile(data_dir, 'BMI_online191119T102404.mat');

strcMask_file   = fullfile(data_dir, 'BMI_roi.mat');

exist(cal_ALL_file)
exist(bmi_file)
exist(strcMask_file)

cal_ALL_data = load(cal_ALL_file); 
bmi_data = load(bmi_file); 
base_data = load(base_file); 
mask_data       = load(strcMask_file)
%%
%Compare how neural data looks in baseline vs BMI: 
%Plot colors: 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 

n_base = ...
    base_data.baseActivity(mask_data.E_base_sel,:);
n_base = ...
    n_base(:, ~isnan(n_base(1,:)));

n_bmi = ...
    bmi_data.data.bmiAct; 
n_bmi = ...
    n_bmi(:, ~isnan(n_bmi(1,:))); 
    
plot_data = n_base; 
t_plot = 1:size(plot_data,2);
[h, offset_vec] = plot_E_activity(t_plot.', plot_data.', mask_data.E_id, E_color, 0);
title('base'); 

plot_data = n_bmi; 
t_plot = 1:size(plot_data,2);
[h, offset_vec] = plot_E_activity(t_plot.', plot_data.', mask_data.E_id, E_color, 0);