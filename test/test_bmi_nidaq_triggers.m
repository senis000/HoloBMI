

%%
s = daq.createSession('ni');
addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
ni_out = [0 0 0]; 
outputSingleScan(s,ni_out);%set   
ni_getimage = [1 0 0]; 
ni_reward   = [0 1 0]; 
ni_holo     = [0 0 1]; 

% %Sound setup:
% xrnd = randn(1000,1);
% reward_sound = audioplayer(xrnd, 10000); %Play sound using: play()
%%

%%
%Sound + Reward
disp('reward delivery!')
play(reward_sound);
pause(0.3); 
outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);

%NOTE: if we do reward sound and solenoid together, solenoid clicks before
%the sound is delivered...

%%
while(true)
    outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);
    pause(0.5)
end

%%
%--------------------------------------------------------------------------
%Do: 
%Check in Prairie View Voltage Recording
%Start wtih Trigger: NONE
%--------------------------------------------------------------------------
disp('test frame trig!')
outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,ni_out);

%%
%--------------------------------------------------------------------------
%Do: 
%Unplug frame trig from AI7, plug in PFI1 (signal for triggering stim)
%Voltage Recording (look at it)
%--------------------------------------------------------------------------
outputSingleScan(s,ni_holo); 
pause(0.01); outputSingleScan(s,ni_out)
disp('done!')
% pause(0.001);outputSingleScan(s,ni_out)

%--------------------------------------------------------------------------
%Do: 
%Replug PFI1, AI7
%Voltage Recording: Start with Trigger: Next Scan
%--------------------------------------------------------------------------
%%
clear s
% pause(0.001);outputSingleScan(s,ni_out)
