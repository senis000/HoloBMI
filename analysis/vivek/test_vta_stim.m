
num_pulses = 28; 
pulse_interval = 71; %ms
pulse_width = 10; 
interpulse = pulse_interval - pulse_width; 

num_samples = 2000; 
stim_tseries = zeros(num_samples,1); 

ind = 1;
for i=1:num_pulses
    stim_tseries(ind:ind+pulse_width) = 1;
    ind = ind + pulse_interval;
end

%%
h = figure;
plot(stim_tseries)

%%
sum(diff(stim_tseries)==-1)