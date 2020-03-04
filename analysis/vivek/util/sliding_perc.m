function [result] = sliding_perc(ts, win, percentile)
%ts: num_obs x num_var
[num_obs num_var] = size(ts); 

num_result_samples = num_obs - (win-1);
result = zeros(num_result_samples, num_var); 

for j = 1:num_var
    for i =1:num_result_samples
        samples_select = i:(i+win-1);
        data_i = ts(samples_select, j);
        result(i, j) = prctile(data_i, percentile);
    end
end
