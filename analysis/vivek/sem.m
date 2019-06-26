function [sem_result, mean_result] = sem(data_over_var_obs, valid_over_var_obs)
%function [sem_result, mean_result] = sem(data_over_var_obs);
%DESCRIPTION:
%Calculates standard error of mean: s/sqrt(n)
%INPUT:
%data_over_var_obs - num_variables X num observations
%valid - num_variables X num_observations

size_check = sum(size(data_over_var_obs) == size(valid_over_var_obs))==2; 
assert(size_check, 'data and valid are not the same size'); 

%%
num_var = size(data_over_var_obs,1);
num_observations = size(data_over_var_obs,2);

mean_result = zeros(num_var, 1);
sem_result = zeros(num_var, 1);
for var_idx = 1:num_var
    valid_idxs  = find(valid_over_var_obs(var_idx, :));
    num_valid   = length(valid_idxs);
    
    data_sel = data_over_var_obs(var_idx, valid_idxs);
    mean_result(var_idx) = ...
        mean(data_sel);
    sem_result(var_idx) = ...
        std(data_sel)/sqrt(num_valid);
end


