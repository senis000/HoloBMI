
%%
addpath(genpath('/Users/vivekathalye/Dropbox/Code/analysis_util')); 
%%

data_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline'; 
exist(data_dir)
baseline_file = fullfile(data_dir, 'baseline_IntegrationRois_00001.csv'); 
exist(baseline_file)
%%
%raw fluorescence 

A = csvread(baseline_file, 1, 0);

%%
ts = A(:,1); 
frame = A(:,2); 
n = A(:, 3:end); 
num_n = size(n,2); 

%
[n_var, I] = sort(var(n, 0, 1), 'descend');
n = n(:,I_nnan); 
n(:,find(isnan(n_var))) = [];
n_var(find(isnan(n_var))) = [];

%%
% 
num_plot = 90;

offset_vec = cumsum(500*ones(num_plot,1)); 

h = figure; 
hold on; 
for i = 1:num_plot
    y_plot = n(:,i)-mean(n(:,i)); 
	plot(ts, y_plot-offset_vec(i)); 
end
% 
% %good channels:
%1, 2, 3, 4, 5, 6, 9, 10

%%
direct = [1 2 3 4 5 6 9 10];
n_direct = n(:,direct); 

E_direct = [1 1 1 1 2 2 2 2]; 
E1 = [1 2 3 4] %goes down, green
E2 = [5 6 9 10] %goes up, blue
num_E1 = length(E1); 
num_E2 = length(E2); 
%%

offset_vec = cumsum(500*ones(num_plot, 1)); 
h = figure;
hold on; 
E_color = {'r', 'k'}
for i=1:length(direct)
    y_plot = n_direct(:,i);
    plot(ts, y_plot-offset_vec(i), E_color{E_direct(i)}); 
end

%%
%Calculate dF/F for each channel: 
n0 = mean(n_direct,1); 
n0_mat = repmat(n0, [size(n_direct,1) 1]); 
size(n0_mat)

dff = (n_direct-n0_mat)./n0_mat;
offset_size = 1.5
offset_vec = cumsum(offset_size*ones(num_plot, 1)); 
h = figure;
hold on; 
E_color = {'r', 'k'}
for i=1:length(direct)
    y_plot = dff(:,i);
    plot(ts, y_plot-offset_vec(i), E_color{E_direct(i)}); 
end

%%
std_dff = var(dff, 0, 1).^(1/2)

dff_z = dff./repmat(std_dff, [size(dff,1) 1]); 
h = figure;
hold on; 
E_color = {'r', 'k'}
offset_size = 10
offset_vec = cumsum(offset_size*ones(num_plot, 1)); 
for i=1:length(direct)
    y_plot = dff_z(:,i);
    plot(ts, y_plot-offset_vec(i), E_color{E_direct(i)}); 
end
title('dffz'); 
xlabel('time'); 
ylabel('dffz'); 

% %%
% h = figure; 
% hold on; 
% 
% plot(ts, dff_z(:,end)); 
% plot(ts, dff_z(:,1)); 

%%
%dff cov: 
dff_z_cov = cov(dff_z);
h = figure;
imagesc(dff_z_cov); 
colorbar
axis square
xlabel('roi')
ylabel('roi')
colormap
caxis([-0.2 0.3]); 

%%
h = figure;
scatter(dff_z(:,3), dff_z(:,4)); 
[rho, pval] = corr(dff_z(:,3), dff_z(:,4));
rho
pval

%%
[u,s,v] = svd(dff_z_cov); 
s_cumsum = cumsum(diag(s))/sum(diag(s)); 
h = figure;
plot(s_cumsum, '.-', 'MarkerSize', 20); 
axis square

% %%
% test_cov = [[ones(4,4) -ones(4,4)]; [-ones(4,4) ones(4,4)]]; 
% h = figure;
% imagesc(test_cov); 
% colorbar
% axis square

% %%
% test_cov = diag(ones(8,1)); 
% h = figure;
% imagesc(test_cov); 
% colorbar
% axis square

%%
decoder = [-ones(num_E1, 1); ones(num_E2,1)]; 
decoder = decoder/norm(decoder); 
cursor_cov = decoder'*dff_z_cov*decoder; 
% cursor_cov = decoder'*test_cov*decoder; 
cursor_cov


%%
% what is the p-value we want? : 
calibrated_rate = 1/(1.5*60*30); %one success per frames in 1.5 minutes
T = norminv((1-calibrated_rate), 0, cursor_cov); 
T

%%
h = figure;
hold on;
scatter(dff_z(:,1), dff_z(:,3));
x = -2:.1:12; 
plot(x, T-x); 

%%
cursor_obs = dff_z*decoder; 
h = figure;
hist(cursor_obs, 20); 
vline(T)

%%
sum(cursor_obs >= T)

%%
T_t = cursor_obs >= T; 
test = dff_z(T_t,:)

h = figure;
imagesc(test)
axis square
colorbar

%%
dff_z_ceil = dff_z; 

target_frac = 0.95; 
for i = 1:size(dff_z, 2)
    ind = (dff_z_ceil(:,i) >= target_frac*T); 
    dff_z_ceil(ind) = target_frac*T;
end

%%
cursor_obs2 = dff_z_ceil*decoder; 
h = figure;
hist(cursor_obs2, 20); 
vline(T)

%%
sum(cursor_obs2 >= T)

%%
%Iterate on T value, until perc correct value is achieved using truncated
%neural activity
frame_rate = 10; 
min_per_reward_range = [1.5 1]; 
reward_per_min_range = 1./min_per_reward_range; 
reward_per_frame_range = reward_per_min_range/(60*frame_rate)

%
cov_bool = 0; 
win_bool = 1;
win_size = 8; %number of samples to average over:
win_filter = ones(win_size,1)/win_size; 

%1) Window the dff:
num_neurons = size(dff_z, 2); 
if(win_bool)
    num_samples = max(length(dff_z)-(win_size-1), 0); 
    %MAX(LENGTH(A)-MAX(0,LENGTH(B)-1),0)
    n_signal = zeros(num_samples, num_neurons); 
    for n_i = 1:num_neurons
        n_signal(:, n_i) = conv(dff_z(:,n_i), win_filter, 'valid'); 
    end
else
    num_samples = size(dff_z,1); 
    n_signal = dff_z; 
end


T0 = max(n_signal*decoder)

%%
T = T0; 
T_delta = 0.01; %0.01*T0; 
T_frac = 0.95;
task_complete = 0;

T_vec = []; 
reward_per_frame_vec = []; 
while(~task_complete)
    T_vec = [T_vec T];
    
    
    %2) Ceiling to T_frac*T:
    n_analyze = n_signal;  
    for n_i = 1:num_neurons
        ind = (n_analyze(:,n_i) >= T_frac*T); 
        n_analyze(ind, n_i) = T_frac*T;
    end
    
    %Calculate p-value of achieving target: 
    if(cov_bool)
        analyze_cov = cov(n_analyze); 
        cursor_cov = decoder'*analyze_cov*decoder; 
        cursor_mean = mean(n_analyze*decoder);
        r_cdf = normcdf(T,cursor_mean, cursor_cov);
        reward_prob_per_frame = 1-r_cdf; %probability of reward per frame
    else
        reward_prob_per_frame = sum((n_analyze*decoder) >= T)/length(n_analyze);
    end
    reward_prob_per_frame
    reward_per_frame_vec = [reward_per_frame_vec reward_prob_per_frame]; 
   
    %Update T:
    if((reward_prob_per_frame >= reward_per_frame_range(1)) && (reward_prob_per_frame <= reward_per_frame_range(2)))
        task_complete = 1;
        disp('task complete!');
    elseif(reward_prob_per_frame > reward_per_frame_range(2))
        T = T+T_delta; 
%         T_delta = 0.1*T;
    elseif(reward_prob_per_frame < reward_per_frame_range(1))
        T = T-T_delta; 
%         T_delta = 0.05*T;        
    end 
end
%Outputs: 
%T
%T_vec
%reward_per_frame
%reward_per_frame_vec
%dff_z_ceil


h = figure;
plot(T_vec); 
xlabel('alg iteration'); 
ylabel('target'); 
if cov_bool
    title('target on data cov'); 
else
    title('target on data'); 
end

h = figure;
hold on; 
plot(reward_per_frame_vec); 
hline(reward_per_frame_range(1)); 
hline(reward_per_frame_range(2)); 
xlabel('alg iteration'); 
ylabel('reward per min'); 
if cov_bool
    title('target on data cov'); 
else
    title('target on data'); 
end

%Plot the actual data: 
cursor_obs = n_analyze*decoder; 
h = figure;
hold on; 
hist(cursor_obs, 50); 
vline(T); 
disp('num baseline hits:'); 
num_baseline_hits = sum(cursor_obs >= T)
if cov_bool
    title('target on data cov'); 
else
    title('target on data'); 
end

%%
h = figure;
scatter(n_analyze(:,1), n_analyze(:,2)); 


%%
%dff vs dff smooth

offset = 5;
offset_vec = 5*(1:num_neurons); 
h = figure; hold on;
for i = 1:num_neurons 
%     plot(n_analyze(:,i)); 
    plot(1:length(dff_z), dff_z(:,i)-offset_vec(i), 'k'); 
    plot(8:length(dff_z), n_signal(:,i)-offset_vec(i), 'r'); 
%     legend({'dff', 'smooth'}); 
end

%%
%smooth vs smooth+ceiling
h = figure; hold on;
for i = 1:num_neurons
%     plot(n_analyze(:,i)); 
    plot(1:length(n_signal), n_signal(:,i)-offset_vec(i), 'k'); 
    plot(1:length(n_analyze), n_analyze(:,i)-offset_vec(i), 'r'); 
    legend({'smooth', 'smooth ceil'});
end
reward_times = find(n_analyze*decoder>=T); 
for i=1:length(reward_times)
    vline(reward_times(i)); 
end