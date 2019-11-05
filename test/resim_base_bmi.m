function [] = resim_base_bmi(data_dir, save_dir, ...
    base_file, ...
    bmi_file, ...    
    roi_data_file, ...
    strcMask_file, ...
    cal_file, ...
    cal_ALL_file)

% 11.4.19
%DEBUG INPUT: 
% roi_data_file   = fullfile(data_dir, 'roi_data.mat'); 
% strcMask_file   = fullfile(data_dir, 'strcMask.mat');
% base_file       = fullfile(data_dir, 'BaselineOnline191103T143329.mat'); 
% bmi_file        = fullfile(data_dir, 'BMI_online191103T162233.mat'); 
% cal_file        = fullfile(data_dir, 'BMI_target_info_20191103T145403.mat');
% cal_ALL_file    = fullfile(data_dir, 'target_calibration_ALL_20191103T145403.mat');  


%Check input paths exist: 
exist(base_file)
exist(bmi_file)
exist(roi_data_file)
exist(strcMask_file)
exist(cal_file)
exist(cal_ALL_file)


mask_data       = load(strcMask_file)
E1_base         = mask_data.E_base_sel(mask_data.E_id==1); 
E2_base         = mask_data.E_base_sel(mask_data.E_id==2);
%
%Constants:
frameRate = 29.989;
