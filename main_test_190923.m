
%main_test_190923
%%
s = daq.createSession('ni');
addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
ni_out = [0 0 0]; 
outputSingleScan(s,ni_out);%set   
ni_getimage = [1 0 0]; 
ni_reward   = [0 1 0]; 
ni_holo     = [0 0 1]; 

%%
%TEST NIDAQ PULSES ARE RECEIVED BY VOLTAGE REC: 
%

%%
%TRIGGER FRAME 
pause_len = 0.01; %0.001
outputSingleScan(s,ni_getimage); pause(pause_len); outputSingleScan(s,ni_out);

%%
%TRIGGER HOLO 
%Put "Monaco Trig" cable into "Monaco Trig REC" to confirm NIDAQ sends
%trigger
%After check, change cables back.  

pause_len = 0.01; %0.001
outputSingleScan(s,ni_holo); pause(pause_len); outputSingleScan(s,ni_out);

%%
% %TRIGGER SOL (Reward / Arduino) Redundant with next test
% outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);

%%
%{
%TRIGGER SOL (Reward / Arduino)    Sound + Reward
0) plug in arduino USB cable and power (needed for sol)
1) Find Ines' sol calibration length for opening solenoid, input the
length into:
'G:\VivekNuria\Code\arduino\sucrose_delivery\sucrose_delivery.ino'.
Verify and Upload (COM12)
2) Run following:
%}
disp('reward delivery!')
% play(reward_sound);

% pause(0.3); 
% for i=1:100
    outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);
    pause(0.3)
% end

%%
%Test audio feedback
task_settings.fb.arduino.com            = 'COM11';
task_settings.fb.arduino.label          = 'Mega2560';
task_settings.fb.arduino.pin            = 'D3';
a = arduino(task_settings.fb.arduino.com, ...
    task_settings.fb.arduino.label);

%%
fb_freq_i = 7000;
task_settings.fb.arduino.duration = 1
playTone(a,...
    task_settings.fb.arduino.pin,...
    fb_freq_i,...
    task_settings.fb.arduino.duration)

%%
clear s
clear a