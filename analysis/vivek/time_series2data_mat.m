function data_mat = time_series2data_mat(time_series, center_times, win)
%7.20.18
%Extract windows of time_series centered at center_times
%time_series: num_var x num_observations 
%win: vector 

% %Debug input:
% time_series = neural;
% center_times = Holo_frames;
% win = [-100 100]

%
num_var = size(time_series, 1); 
num_samples = win(2)-win(1)+1;
num_trials = length(center_times);

data_mat = zeros(num_var, num_samples, num_trials); 
for trial_i=1:length(center_times)
    time_i = center_times(trial_i);
    win_i = time_i + win;
    data_mat(:,:,trial_i) = time_series(:,win_i(1):win_i(2));
end