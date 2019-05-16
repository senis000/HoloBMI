

%%
s = daq.createSession('ni');
addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
ni_out = [0 0 0]; 
outputSingleScan(s,ni_out);%set   
ni_getimage = [1 0 0]; 
ni_reward   = [0 1 0]; 
ni_holo     = [0 0 1]; 

%Sound setup:
% xrnd = randn(1000,1);
% reward_sound = audioplayer(xrnd, 10000); %Play sound using: play()
%%
%Sound + Reward
% weight of the black cup for reward calibration is: 1.596

disp('reward delivery!')
t = [];
for i=1:100
%     tic;play(reward_sound);t(end+1) = toc;
%     tic;sound(xrnd, 10000);t(end+1) = toc;
    pause(0.5);
    outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);
end
disp ('end')
% pause(0.3); 
% outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);

%NOTE: if we do reward sound and solenoid together, solenoid clicks before
%the sound is delivered...

%%
disp('test frame trig!')
outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,ni_out);



%%
outputSingleScan(s,ni_holo); 
pause(0.01); outputSingleScan(s,ni_out)
disp('done!')
% pause(0.001);outputSingleScan(s,ni_out)

