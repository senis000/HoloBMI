function [fb_settings] = define_fb_audio_settings()

fb_settings.target_low_freq        = 1; 
%Set the target cursor value to be the low frequency
fb_settings.freq_min               = 6000; 
fb_settings.freq_max               = 19000; 
fb_settings.arduino.com            = 'COM11';
fb_settings.arduino.label          = 'Mega2560';
fb_settings.arduino.pin            = 'D3';
fb_settings.arduino.duration       = 0.3; %ms, tones update at rate of BMI code, this is the longest a tone will play for
fb_settings.min_perctile            = 10; 
