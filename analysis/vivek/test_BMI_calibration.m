%test_BMI_calibration
clear all
%%
addpath(genpath('/Users/vivekathalye/Dropbox/Code/analysis_util')); 
%%
data_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline'; 
exist(data_dir)
% baseline_file = fullfile(data_dir, 'baseline_IntegrationRois_00001.csv'); 
baseline_file = fullfile(data_dir, 'baseline_IntegrationRois_00001_1.csv'); 
exist(baseline_file)

%%
%raw fluorescence 
A = csvread(baseline_file, 1, 0);
ts = A(:,1); 
frame = A(:,2); 
n = A(:, 3:end); 
num_n = size(n,2); 

[n_var, I] = sort(var(n, 0, 1), 'descend');
n = n(:,I); 
n(:,find(isnan(n_var))) = [];
n_var(find(isnan(n_var))) = [];

num_plot = 20;

offset = 200; 
offset_vec = cumsum(offset*ones(num_plot,1)); 

h = figure; 
hold on; 
for i = 1:num_plot
    plot_i = i+20
    y_plot = n(:,plot_i)-mean(n(:,plot_i)); 
	plot(ts, y_plot-offset_vec(i)); 
end

%% 
% baseline_IntegrationRois_00001_1.csv
direct = [1 2 4 5 6 17 27 39]
n_direct = n(:,direct); 
E_direct = [2 2 2 2 1 1 1 1]'; 

% %--------------------------------------------------------------------------
% % baseline_IntegrationRois_00001.csv
% direct = [1 2 3 4 5 6 9 10];
% n_direct = n(:,direct); 
% E_direct = [1 1 1 1 2 2 2 2]'; 

%%
E1_sel = E_direct==1; 
E1_sel_idxs = find(E1_sel); 
E2_sel = E_direct==2; 
E2_sel_idxs = find(E2_sel); 
E1 = direct(E1_sel); %goes down
E2 = direct(E2_sel); %goes up
num_E1 = length(E1); 
num_E2 = length(E2); 
num_neurons = num_E1 + num_E2; 

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
offset_size = 4
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
offset_size = 15
offset_vec = cumsum(offset_size*ones(num_plot, 1)); 
for i=1:length(direct)
    y_plot = dff_z(:,i);
    plot(ts, y_plot-offset_vec(i), E_color{E_direct(i)}); 
end
title('dffz'); 
xlabel('time'); 
ylabel('dffz'); 

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
xlabel('PC'); 
ylabel('Frac Var Explained'); 

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

E1_proj = zeros(num_neurons, 1); 
E1_proj(E1_sel) = 1;
E1_norm = sum(E1_sel); %can replace with vector norm

disp('E1 proj'); 
E1_proj = E1_proj/E1_norm;

E2_proj = zeros(num_neurons, 1); 
E2_proj(E2_sel) = 1; 
E2_norm = sum(E2_sel); 
disp('E2 proj')
E2_proj = E2_proj/E2_norm

disp('decoder:')
decoder = E2_proj - E1_proj
% decoder = decoder/norm(decoder); 
cursor_cov = decoder'*dff_z_cov*decoder; 
% cursor_cov = decoder'*test_cov*decoder; 
cursor_cov

%%
E1_signal = dff_z*E1_proj;
E2_signal = dff_z*E2_proj;
h = figure;
scatter(E1_signal, E2_signal);
axis square
xlabel('E1'); 
ylabel('E2'); 

%%
%Calibrate the BMI: 
% 1. E2-E1>a: Linear constraint is met
% 2. nontop_E2 > b: One E2 neuron doesn?t do the work. The percentile of the non-top neurons is above a threshold, so the net effect of non-top neurons is still positive
% 3. E1<c: E1 does not go up too much.  
% Three parameters: a, b, c
% a = fit to chance rate on {covariance / real data}
% b = 50th percentile
% c = 50th percentile


%%
%Iterate on T value, until perc correct value is achieved using truncated
%neural activity
frame_rate = 10; 
min_per_reward_range = [1.5 1]; 
reward_per_min_range = 1./min_per_reward_range; 
reward_per_frame_range = reward_per_min_range/(60*frame_rate);

%
cov_bool = 0; 
win_bool = 1;
win_size = 2; %number of samples to average over:
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

%
%Real data
n_analyze = n_signal;
analyze_cov = cov(n_analyze);
analyze_mean = mean(n_analyze); 
E1_mean = mean(analyze_mean(E1_sel));
E2_subord_mean = zeros(num_E2,1);
E2_subord_std = zeros(num_E2,1); 
E1_analyze = n_analyze(:,E1_sel); 
E2_analyze = n_analyze(:,E2_sel); 
for E2_i = 1:num_E2
    subord_sel = E2_sel;
    subord_sel(E2_sel_idxs(E2_i)) = 0; 
    
    E2_subord_mean(E2_i) = mean(analyze_mean(subord_sel));     
%     E2_subord_mean(E2_i) = sum(analyze_mean(subord_sel))-analyze_mean(E2_sel_idxs(E2_i));
%     %renormalize:
%     E2_subord_mean(E2_i) = E2_subord_mean(E2_i)/(num_E2-1);

    var_i = subord_sel'*analyze_cov*subord_sel; 
    E2_subord_std(E2_i) = sqrt(var_i);     
end

E2_sum_analyze = sum(E2_analyze,2); 

%signals needed for target detection:
cursor_obs                      = n_analyze*decoder; 
E1_mean_analyze                 = mean(E1_analyze,2);
[E2_dom_samples, E2_dom_sel]    = max(E2_analyze, [], 2);
E2_subord_mean_analyze          = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);


%%
T0 = max(n_analyze*decoder);
min_cursor = min(n_analyze*decoder);

E2_coeff = 0.5; 
E2_subord_thresh = E2_subord_mean+E2_subord_std*E2_coeff;

T = T0; 
T_delta = 0.05; 
E2_coeff_delta = 0.05; 
task_complete = 0;
T_vec = []; 
reward_per_frame_vec = []; 
 
rand_num_samples = 500000;
while(~task_complete)
    T_vec = [T_vec T];
    if(cov_bool)
        %Generate samples
        mvn_samples = mvnrnd(analyze_mean,analyze_cov,rand_num_samples);
        
        %Check probability of success condition
        cursor_samples = mvn_samples*decoder; 
        E1_samples = mvn_samples(:,E1_sel); 
        E2_samples = mvn_samples(:,E2_sel); 
        E1_mean_samples = mean(E1_samples,2); 
        E2_mean_samples = mean(E2_samples,2); 
        
        %1) E2-E1 > alpha
        c1 = find(cursor_samples>= T); 
        
        %2) E1 < mu
        c2 = find(E1_mean_samples <= E1_mean); 
        
        %3) E2_subord > mu (anded with previous constraint)
        [E2_dom_samples, E2_dom_sel_samples] = max(E2_samples, [], 2);
        E2_subord_sum_samples = E2_mean_samples*num_E2 - E2_dom_samples;
        E2_subord_mean_samples = E2_subord_sum_samples/(num_E2-1); 
        %For each idx, subtract the 
        c3 = find(E2_subord_mean_samples >= E2_subord_thresh(E2_dom_sel_samples)); 
        %E2_subord_mean(E2_dom_sel_samples)); 
        
        %Intersect the indices:
        hits = intersect(intersect(c1, c2), c3); 
        reward_prob_per_frame = length(hits)/rand_num_samples; 
    else
        %1) E2-E1 > alpha
        c1 = find(cursor_obs >= T); 
        %2) E1 < mu
        c2 = find(E1_mean_analyze <= E1_mean);
        %3) E2_subord > mu (anded with previous constraint)
        %For each idx, subtract the 
        c3 = find(E2_subord_mean_analyze >= E2_subord_thresh(E2_dom_sel)); 

        hits = intersect(intersect(c1, c2), c3);        
        reward_prob_per_frame = length(hits)/length(n_analyze);
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
        if(T<=0)
            T=T0; 
            E2_coeff = E2_coeff - E2_coeff_delta;
            E2_subord_thresh = E2_subord_mean+E2_subord_std*E2_coeff;
%             T_delta = 0.5*T_delta;     
        end        
        T = T-T_delta; 
    end
    T
end 

%%
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


%%
%Plot the actual data: 
disp('T'); 
T
%data: 2.29
cursor_obs = n_analyze*decoder; 
c1 = find(cursor_obs >= T); 
disp('num E2-E1 >= T'); 
length(c1)

E1_mean_analyze = mean(E1_analyze,2)
c2 = find(E1_mean_analyze <= E1_mean);
disp('E1 >= b'); 
length(c2)

[E2_dom_samples, E2_dom_sel] = max(E2_analyze, [], 2);
E2_subord_mean_analyze = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);
%For each idx, subtract the 
c3 = find(E2_subord_mean_analyze >= E2_subord_thresh(E2_dom_sel)); 
disp('E2 subord >= c'); 
length(c3)

hit_times = intersect(intersect(c1, c2), c3);
disp('num baseline hits:'); 
num_hits = length(hit_times)

%%
cursor_obs = n_analyze*decoder; 
h = figure;
hold on; 
hist(cursor_obs, 50); 
vline(T); 
disp('num baseline hits:'); 

if cov_bool
    title('target from cov'); 
else
    title('target from data'); 
end

%%
%Plot the hit times: 
h = figure; hold on;
E_color = {'r', 'k'}
offset = 10;
offset_vec = cumsum(offset*ones(length(E_direct), 1)); 
for i = 1:num_neurons
    plot(1:length(n_analyze), n_analyze(:,i)-offset_vec(i), E_color{E_direct(i)}); 
end

%c1:
c1_offset = offset_vec(end)+offset;
plot(1:length(cursor_obs), cursor_obs-c1_offset);
hline(T-c1_offset)

%c2:
c2_offset = offset_vec(end)+2*offset;
plot(1:length(E1_mean_analyze), E1_mean_analyze-c2_offset);
hline(E1_mean-c2_offset)

%c3:
c3_offset = offset_vec(end)+3*offset;
plot(E2_subord_mean_analyze-c3_offset); 
plot(E2_subord_thresh(E2_dom_sel)-c3_offset);

for i=1:length(hit_times)
    vline(hit_times(i)); 
end

%%
h = figure;
plot(E2_subord_mean(E2_dom_sel));


%%
cursor_obs = n_analyze*decoder; 
c1 = find(cursor_obs >= T); 
disp('num E2-E1 >= T'); 
length(c1)

E1_mean_analyze = mean(E1_analyze,2)
c2 = find(E1_mean_analyze <= E1_mean);
disp('E1 >= b'); 
length(c2)

[E2_dom_samples, E2_dom_sel] = max(E2_analyze, [], 2);
E2_subord_mean_analyze = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);
%For each idx, subtract the 
c3 = find(E2_subord_mean_analyze >= E2_subord_mean(E2_dom_sel)); 
disp('E2 subord >= c'); 
length(c3)

%%
h = figure;
hold on;
plot(1:length(cursor_obs), cursor_obs); 
hline(T); 
