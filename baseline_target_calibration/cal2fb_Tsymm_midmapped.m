function [cal] = cal2fb_Tsymm_midmapped(cal, task_settings)
%INPUT: 
%cal
%FIELDS: 
% -- target.T
% -- fb.freq_min
% -- fb.fre_max

%Copy from task_settings:
cal.fb.fb_bool = ...
    task_settings.fb.fb_bool; 
cal.fb.target_low_freq = ...
    task_settings.fb.target_low_freq; 
cal.fb.freq_min = ...
    task_settings.fb.freq_min; 
cal.fb.freq_max = ...
    task_settings.fb.freq_max; 

%Calculate:
cal.fb.cursor_min         = ...
    cal.target.E1_hit_cal.T; 
%negate because E1 was calibrated to E1mE2
cal.fb.cursor_max         = ...
    cal.target.E2_hit_cal.T; %for fb, cursor is floor to this value
cal.fb.cursor_range       = ...
    cal.fb.cursor_max - cal.fb.cursor_min; 
% freq = a*exp(b*(cursor_trunc-cursor_min))

cal.fb.a                  = ...
    cal.fb.freq_min; 
cal.fb.b                  = ...
    (log(cal.fb.freq_max) - log(cal.fb.a))/cal.fb.cursor_range; 

%Divide into Lower Half, Upper Half mapping of cursor to frequency
cal.fb.cursor_middle    = ...
    cal.target.cursor_middle;

cal.fb.freq_middle    = ...
    cal.fb.a*exp(cal.fb.b*cal.fb.cursor_range/2);    
%cursor2audio_freq(cal.fb.cursor_middle, cal);

cal.fb.a_bin = zeros(1,2); 
cal.fb.b_bin = zeros(1,2); 


if cal.fb.target_low_freq == 1
    %This means cursor up makes auditory freq go down:
    cursor_min      = -cal.fb.cursor_max;
    cursor_middle   = -cal.fb.cursor_middle; 
    cursor_max      = -cal.fb.cursor_min;
    
else
    %This means cursor up makes auditory freq go up:
    cursor_min      = cal.fb.cursor_min;
    cursor_middle   = cal.fb.cursor_middle; 
    cursor_max      = cal.fb.cursor_max;    
end

%Lower Half:
cal.fb.a_bin(1) = ...
    cal.fb.freq_min;
cal.fb.b_bin(1) = ...
    (log(cal.fb.freq_middle) - log(cal.fb.a_bin(1)))/(cursor_middle - cursor_min); 
%Upper Half:
cal.fb.a_bin(2) = ...
    cal.fb.freq_middle;
cal.fb.b_bin(2) = ...
    (log(cal.fb.freq_max) - log(cal.fb.a_bin(2)))/(cursor_max - cursor_middle); 

end