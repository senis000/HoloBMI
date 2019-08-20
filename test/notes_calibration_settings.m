calibration_settings

TODO: organize baseline calibration to save the following params
TODO: write function for auditory feedback

%PATHS
BMI_roi_path
plotPath

%REWARD 
%calculate with baseline file and task_settings telling you baseline_len
num_base_samples          = []; %number of samples collected in baseline
baseline_frameRate        = []; 
frames_per_reward_range   = []; %converts sec_per_reward_range

%BMI neurons
strcMask
E1_base
E2_base
E_base_sel
E_id
E1_sel_idxs
E2_sel_idxs

%Decoder
decoder
E1_proj
E1_norm
E2_proj
E2_norm

%Target
T
E1_thresh
E2_subord_thresh
%Target related
E1_mean
E1_std
E2_subord_mean 
E2_subord_std
n_mean
n_std

%Cal results
num_c1
num_c2
num_c3
num_hits_no_b2base
num_valid_hits

%Auditory feedback



%TODO: 
%in baseline, use 'task_settings'

hits = intersect(intersect(c1, c2), c3);

hit_times_no_b2base = intersect(intersect(c1, c2), c3);