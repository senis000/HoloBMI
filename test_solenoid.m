

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