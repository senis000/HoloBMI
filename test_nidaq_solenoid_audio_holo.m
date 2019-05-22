

%%
s = daq.createSession('ni');
addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
ni_out = [0 0 0]; 
outputSingleScan(s,ni_out);%set   
ni_getimage = [1 0 0]; 
ni_reward   = [0 1 0]; 
ni_holo     = [0 0 1]; 
    
%%
disp('reward delivery!')
outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);

%%
%To Observe this working: 
%1) Open Prairie
%2) Sync: Voltage Recording (Current), Mark Points (Current)
%3) Mark Points: 
%Add a single point to the point series
%Put laser power >0
%Set iterations > 1
%Wait for Trigger: First Repitition
%Trigger Selection: Start with PFI1

outputSingleScan(s,ni_holo); 
pause(0.01); outputSingleScan(s,ni_out)
disp('done!')
% pause(0.001);outputSingleScan(s,ni_out)

%%
disp('test frame trig!')
outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,ni_out);

% %%
% xrnd = randn(1000,1);
% reward_sound = audioplayer(xrnd, 10000); %Play sound using: play()
% 
% %%
% play(reward_sound)

