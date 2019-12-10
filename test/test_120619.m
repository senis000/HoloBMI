%%
%Goal: 
%choose subset of feedback signal to playback to the animal.
%must meet the following criteria: 
%number of rewards in the interval on both sides must be equal
%the cursor distribution should be symmetric

%A good one for tone symmetry but the rewards (10 vs 12)
%load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-11-16\BMI_cal_ALL_20191116T091949.mat')

%
%%
load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-11-16\BMI_cal_ALL_20191116T091949.mat')

%%
load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-11-22\BMI_cal_ALL_20191122T112510.mat'); 

%%
%Test load just the cursor observations: 

base_cursor_data = ...
    load('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-11-16\BMI_cal_ALL_20191116T091949.mat', 'cursor_obs').cursor_obs;
%load(cal.paths.cal_all, 'cursor_obs').cursor_obs; 

%%
%Test load the BMI online data: 
bmi = ...
    load(fullfile('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-11-16', 'BMI_online191116T100508.mat'))

%%
save_dir = 'D:\Dropbox\Data\bmi_striatum\playback_baseline'

%%

h = figure;

E2_T                = cal.target.E2_hit_cal.T;
E2_hits_b2base      = cal.target.E2_hit_cal.num_hits_b2base;
E2_hits_no_b2base   = cal.target.E2_hit_cal.num_hits_no_b2base;

E1_T                = cal.target.E1_hit_cal.T;
E1_hits_b2base      = cal.target.E1_hit_cal.num_hits_b2base;
E1_hits_no_b2base   = cal.target.E1_hit_cal.num_hits_no_b2base;

h = figure;
hold on; 
hist(cursor_obs, 50); 
vline(E2_T); 
vline(E1_T); 

xlabel('Cursor'); 
ylabel('Number of Observations'); 

title(['E2: T: ' num2str(E2_T) ' hits b2base: ' num2str(E2_hits_b2base) ' hits no b2b: ' num2str(E2_hits_no_b2base) ...
    ' E1: T: ' num2str(E1_T) ' hits b2base: ' num2str(E1_hits_b2base) ' hits no b2b: ' num2str(E1_hits_no_b2base)]); 
% 
%%
fb_obs = cursor2audio_freq_middle_match(cursor_obs, cal)
num_fb_bins = 100; 
h = figure;
hist(fb_obs, num_fb_bins); 
xlabel('audio freq'); 
ylabel('baseline counts'); 
% saveas(h, fullfile(plotPath, 'base_freq_hist.png')); 

%%
fb_obs = cursor2audio_freq_middle_match(cursor_obs, cal)
num_fb_bins = 100; 
h = figure;
hist(log(fb_obs), num_fb_bins); 
xlabel('audio freq'); 
ylabel('baseline counts'); 

%%
%Choose an X sample segment
sample_len_time = 7.5
frame_period = 1000/30 %33 ms
sample_len_time*60*frame_period
num_samples_sel = round(sample_len_time*60*frame_period)

%
sample_offset_start       = 0; 
sel_idxs_start            = (1:num_samples_sel) + sample_offset_start; 
c_sel_start               = cursor_obs(sel_idxs_start); 

sample_offset_end           = length(cursor_obs)-num_samples_sel; 
sel_idxs_end                = (1:num_samples_sel) + sample_offset_end; 
c_sel_end                   = cursor_obs(sel_idxs_end);

%%
sum(c_sel_start >= E2_T)
sum(c_sel_start <= E1_T)

%%
sum(c_sel_end >= E2_T)
sum(c_sel_end <= E1_T)

%%
N_start = hist(c_sel_start, 50); 
h = figure;
hist(c_sel_start, 50); 
vline(E2_T); 
vline(E1_T); 
title('cursor sel start'); 

N_end = hist(c_sel_start, 50); 
h = figure;
hist(c_sel_end, 50); 
vline(E2_T); 
vline(E1_T); 
title('cursor sel end'); 
%%
%We want the segment to be symmetric and have roughly the same number of
%rewards on both sides.  

%TODO: simulate the number of rewards
%For now, just plot the histogram: 


c_sel           = c_sel_start;  

h = figure;
hold on; 
% hist(cursor_obs, 50); 
hist(c_sel, 50); 

vline(E2_T); 
vline(E1_T); 

xlabel('Cursor'); 
ylabel('Number of Observations'); 

title(['E2: T: ' num2str(E2_T) ' E1: T: ' num2str(E1_T) ]); 

%%
%Save out the results: 
%1) cursor_obs
%2) fb_obs
%3) cal

save_path = fullfile(save_dir, 'playback_base111619.mat'); 
save(save_path, 'cursor_obs', 'fb_obs', 'cal'); 

%%
test = load(save_path)

%%
tic
E1_bool = 0
E_cal = cal.target.E2_hit_cal
[E_hit_cal, E_hit_data] = compute_hits_b2base(E1_bool, E_cal, ...
    cal.target.T_prctile, cursor_obs, ...
    [], [], ...
    task_settings.b2base_coeff, task_settings.back2BaseFrameThresh);
toc
E_hit_cal

%%
tic
    E1_T    = prctile(cursor_obs, 100-cal.target.T_prctile); 
    E2_T    = prctile(cursor_obs, cal.target.T_prctile); 
    mid_T   = prctile(cursor_obs, task_settings.fb.middle_prctile); 
toc
%%


[E_hit_cal, E_hit_data] = compute_hits_b2base(E1_bool, E_cal, ...
    T_prctile, cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh)

%%
%need function to do simulation of BMI rapidly to calculate number of
%rewards