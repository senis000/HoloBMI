function [task_settings] = define_BMI_task_settings(fb_bool)

task_settings.frameRate = 29.989; 
%Imaging environment file for baseline acquisition
task_settings.baseline_env = ...
    fullfile('G:\VivekNuria\utils', "Tseries_VivekNuria_15.env");
%Imaging environment file for BMI acquisition
task_settings.bmi_env = ...
    fullfile('G:\VivekNuria\utils', "Tseries_VivekNuria_40.env");
%TODO: BaselineAcqnvsPrairie to take this as input
%TODO: BMIAcqnvsPrairienoTrialsHoloCL_fb_debug_enable to take this as input


%calibration: 
task_settings.calibration.target_on_cov_bool     = 0; 
task_settings.calibration.sec_per_reward_range   = [120 80]; 
task_settings.calibration.baseline_len           = 15*60; %seconds
task_settings.calibration.f0_win_bool            = 1; %during calibration, 
%estimate f0 using the window
task_settings.calibration.f0_init_slide          = 0; %during calibration, 
%if 0, f0 is only used after f0_win samples.  if 1, f0 is
%adapted in the window from 0 to f0_win samples.
%Fields that are set in calibration: 
% task_settings.calibration.num_base_samples          = []; %number of samples collected in baseline
% task_settings.calibration.baseline_frameRate        = []; 
% task_settings.calibration.frames_per_reward_range   = []; %converts sec_per_reward_range

%bmi: 
task_settings.onacid_bool               = false; 
task_settings.bmi_len                   = 40*60; %seconds
task_settings.prefix_win                = 100; 
%set this to 'nonBufferUpdateCounter', 'initFrameBase', number samples to ignore at start of bmi acqn

task_settings.f0_win                    = 2*60*ceil(task_settings.frameRate); 
%'baseFrames' in BMI code, number samples to use to estimate f0

task_settings.dff_win                   = 4; 
%'movingAverageFrames', number of frames to use for smoothing dff

task_settings.cursor_zscore_bool        = 0; 
%- if 1, neural activity is zscored before going into
%cursor calculation. if 0, neural activity is not zscored.  

task_settings.rewardDelayFrames         = 10; 
%TODO confirm arduino code triggers reward immediately, so it doesn't add
%extra delay

task_settings.back2BaseFrameThresh      = 2; %2 frames of back2base 
task_settings.relaxationTime            = 0; 
task_settings.b2base_coeff              = 0.5; 
%a frame counts as back2base if cursor < b2base_coeff*T, where T is target cursor value.  

%feedback: 
task_settings.fb.fb_bool                = fb_bool; 
task_settings.fb.target_low_freq        = 1; 
%Set the target cursor value to be the low frequency
task_settings.fb.freq_min               = 6000; 
task_settings.fb.freq_max               = 19000; 
task_settings.fb.arduino.com            = 'COM11';
task_settings.fb.arduino.label          = 'Mega2560';
task_settings.fb.arduino.pin            = 'D3';
task_settings.fb.arduino.duration       = 0.3; %ms, tones update at rate of BMI code, this is the longest a tone will play for

%Fields set in calibration: 
% task_settings.fb.cursor_min         = []; %for fb, cursor is ceil to this value
% task_settings.fb.cursor_max         = []; %for fb, cursor is floor to this value
% task_settings.fb.cursor_range       = []; 
% % freq = a*exp(b*(cursor_trunc-cursor_min))
% task_settings.fb.a                  = []; 
% task_settings.fb.b                  = []; 

%%
%ToDo data to baseline calibration, BMI acqn: 
%roi_data 
%roi_ctr

%get rid of 'dff_win_bool' from calibration

%%
%Baseline inputs: 
%baseline_frameRate

%%
%Baseline outputs: 
%frames_per_reward_range
%

