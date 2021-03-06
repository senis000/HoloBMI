%test_clda_121319

%%
%{
Notes from mainProt: 
In future: zero z once we've chosen an imaging plane. (FOV)
TODO: debug_input for BMI_CLDA (to test it)
TODO: Change baseline code to be able to use the CLDA data to do a new baseline
TODO: check if deltaf/f range changes over teh course of training.
TODO: why did calibrating on full clda data lead me to a harder threshold
than what was calculated online?  Why is plot_CLDA giving me lower number
of hits than full recalibration? 
TODO: update implementation of F0 to be a percentile of past window.
TODO: code to automate the neuron selection for BMI

TODO: implement init F0 for BMI_CLDA
%}

%%
%{
mainProt code sequence: 
1 baseline with playback
2 baseline calibration 
3 clda (online percentile evaluation)
4 plot online percetile evaluation
5 compare online percentile evaluation to full recalibration

1
baseline with playback: 
    BaselineAcqnvsPrairie_fb_playback
    [baseline_mat, baseline_dat] = ...
        BaselineAcqnvsPrairie_fb_playback(folder, animal, day, AComp, roi_mask, task_settings, ...
        a, playback_data, playback_bool);

2
baseline calibration on the playback
    % Manually adjust task parameters: 
    task_settings.calibration.baseline_len              = 7.5*60; %seconds
    task_settings.calibration.sec_per_reward_range      = [50 35]
    [cal, BMI_roi_path] = baseline2two_target_linear_fb_no_constraint(n_f_file, roi_data_file, task_settings, ...
        E1_base, E2_base, savePath);

3
clda block, evaluates the percentiles (chosen from calibration to playback)
with clda data, in closed loop.
    [clda_mat] = BMI_CLDA(folder, animal, day, ...
        expt_str, cal, task_settings, a, vectorHolo, vectorVTA, ...
        base_cursor_samples, ...
        debug_bool, debug_input);

4
plot the clda results 
[cal_update] = plot_CLDA(clda_mat, clda_dir);

5
compare online percentile evaluation to batch re-calibration on the
clda data
    [cal, BMI_roi_path] = baseline2two_target_linear_fb_no_constraint_input_data(...
        clda_base, ...
        roi_data_file, task_settings, ...
        E1_base, E2_base, savePath);


Reminder:
baseline cal happens with tone playback

To verify: 
-is the clda code updating the percentiles correctly online?  
--need a testbench for the clda code.  
-how many hits does the final clda_cal achieve if it were applied to the clda
block? 
-how many hits would 'cal_init' achieve in clda block? 

Decisions: 
does it make sense to use the short random playback to choose a percentile?
that percentile is determined A) over a short time period, and B) not in
closed loop, so it is unlikely to be a good target percentile for a long
BMI session.  
Possible benefit: if for some reason the neurons need an absurd percentile
(either too low or too high)
to achieve the target reward range, then those neurons may not be good
choices for BMI.  
So the percentile needed to achieve reward chance rate is an indicator of
good neurons to choose.  

An alternative is that we choose a fixed percentile to set as the target just
for the clda.  
Argument: it doesn't matter too much what the percentile is anyway.  
We will choose the exact percentile after the clda block

%}

%%
%{
Code Notes

We need a fast function that simulates BMI given different types of inputs:
1) F, cal
2) F0, F, cal
3) dff, cal
4) cursor, cal
...


Calculates cursor and hits with each dff sample.  
[cursor, E2_hit, E1_hit] = ...
    dff2cursor_target_Esymm_no_constraint(dff, cal)
%NOTE: 
%dff must be of size 1xnum_neurons, or num_neuronsx1


%}

%%
%Files: 
% data_dir = 'Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-13'
data_dir = '/Users/vivekathalye/Dropbox/Data/bmi_striatum/clda_test_121319'
base_cal = load(fullfile(data_dir, 'BMI_cal_ALL_20191213T122651.mat')); 

%%
%{
%NOTE: total recalibration involves re-normalizing data by its range, which
%can change a lot...
can we do analysis to see if the range of a neuron's deltaF/F changes ? 

%}
%CLDA: 
clda_cal = load(fullfile(data_dir, 'CLDA_191213T112010.mat')); 

%%
%Baseline: 

%%
cal_init.target.E2_hit_cal

%%
valid_idxs      = find(~isnan(data.E2_T)); 
last_valid_idxs = valid_idxs(end)
%%
data = clda_cal.data; 
E2_T_valid      = data.E2_T(valid_idxs); 
E1_T_valid      = data.E1_T(valid_idxs); 
mid_T_valid     = data.mid_T(valid_idxs); 
h = figure; 
plot(E2_T_valid); 
title('E2 T'); 

h = figure;
plot(E1_T_valid); 
title('E1 T'); 

% %%
% h = figure; hold on;
% plot(E2_T_valid);
% plot(E1_T_valid); 
% % h = figure;
% % plot(mid_T_valid); 
% % title('Mid T');

%%
cursor_valid_idxs = find(~isnan(clda_cal.data.cursor)); 
last_cursor_valid_idx = cursor_valid_idxs(end)

%%
data_test = clda_cal.

