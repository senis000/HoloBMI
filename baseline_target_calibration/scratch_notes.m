% input params to baseline2target_decoder_vE1strict

n_f_file, Acomp_file, onacid_bool,  ...
    E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
    prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir, ...
    cursor_zscore_bool, f0_init_slide

--
auditory_fb_bool 
freq_min
freq_max
freq_target (assert it is either freq_min, freq_max)
cursor_min_perc (which cursor percentile is mapped to the min freq)


--
freq_target:
(assert it should be either min or max)
if it's min, negate cursor.  
    
--
Inputs to baseline function: 

neural
--
n_f_file
Acomp_file
onacid_bool
E1_base
E2_base

calibration
--
frames_per_reward_range
target_on_cov_bool

bmi
--
prefix_win
f0_win_bool
f0_win
dff_win_bool
dff_win
cursor_zscore_bool
f0_init_slide

cursor2freq
--
auditory_fb_bool 
freq_min
freq_max
freq_target (assert it is either freq_min, freq_max)
cursor_min_perc (which cursor percentile is mapped to the min freq)

save_dir















--
manual drawing of spatial mask for each neurons


