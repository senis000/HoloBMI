function [cal_update] = plot_BMI_online(bmi_online_mat)

%%
%Plot the CLDA parameter evolution
%
%Plot histogram with final settings
%
%
%%
load(bmi_online_mat); 

%%
%Get the correct number of hits with the final calibration parameters: 
cal_update = cal; 

E1_bool     = 1; 
E_cal       = cal.target.E1_hit_cal; 
T_prctile   = 100-cal.target.T_prctile;
cursor_obs  = data.cursor; 
mean_E2     = []; 
mean_E1     = []; 

b2base_coeff        = task_settings.b2base_coeff;
b2baseFrameThresh   = task_settings.back2BaseFrameThresh;

[E1_hit_cal, E1_hit_data] = compute_hits_b2base(E1_bool, E_cal, ...
    T_prctile, cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh);
cal_update.E1_hit_cal = E1_hit_cal; 

%%
E1_bool     = 0; 
E_cal       = cal.target.E2_hit_cal; 
T_prctile   = cal.target.T_prctile;
cursor_obs  = data.cursor; 
mean_E2     = []; 
mean_E1     = []; 

b2base_coeff        = task_settings.b2base_coeff;
b2baseFrameThresh   = task_settings.back2BaseFrameThresh;

[E2_hit_cal, E2_hit_data] = compute_hits_b2base(E1_bool, E_cal, ...
    T_prctile, cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh);
cal_update.E2_hit_cal = E2_hit_cal; 

%%
%Summary results: 
disp('E2mE1 T:'); 
cal_update.E2_hit_cal.T

disp('E2 HIGH: num valid hits (WITH B2BASE):'); 
cal_update.E2_hit_cal.num_hits_b2base

disp('E2 HIGH: num baseline hits WITHOUT B2BASE:'); 
cal_update.E2_hit_cal.num_hits_no_b2base

disp('E1 HIGH: num valid hits (WITH B2BASE):'); 
cal_update.E1_hit_cal.num_hits_b2base

disp('E1 HIGH: num baseline hits WITHOUT B2BASE:'); 
cal_update.E1_hit_cal.num_hits_no_b2base

%%

