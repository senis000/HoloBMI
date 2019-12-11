%test_clda_121119

clda = load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-10\CLDA_191210T120802.mat')

%%
%
bmi = load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-10\BMI_online191210T122305.mat'); 

%%
base = load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-10\BMI_cal_ALL_20191210T120615.mat')

%%
valid_idxs      = ~isnan(data.E2_T); 
E2_T_valid      = data.E2_T(valid_idxs); 
E1_T_valid      = data.E1_T(valid_idxs); 
mid_T_valid     = data.mid_T(valid_idxs); 
h = figure;
plot(E2_T_valid); 
title('E2 T'); 

h = figure;
plot(E1_T_valid); 
title('E1 T'); 

h = figure;
plot(mid_T_valid); 
title('Mid T'); 

%%
disp('init E1 cal'); 
cal_init.target.E1_hit_cal

disp('init E2 cal'); 
cal_init.target.E2_hit_cal

%%
disp('init E1 cal'); 
cal.target.E1_hit_cal

disp('init E2 cal'); 
cal.target.E2_hit_cal

%%
bmi.cal.target.E2_hit_cal

%%
cursor_valid = data.cursor; 
cursor_valid(isnan(cursor_valid)) = []; 
h = figure;
hist(cursor_valid, 50); 
vline([cal.target.E2_hit_cal.T E2_T_valid(end)]); 

%%
prctile(cursor_valid, cal.target.T_prctile)

%%
prctile(base.cursor_obs, cal.target.T_prctile)