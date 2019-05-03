function baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
%4.18.19
%inputs:
%n_f_file - contains matrix, neural fluorescence from baseline file, num_samples X num_neurons_baseline 
%Acomp_file - 
%E1_base - E1 idxs in baseline data
%E2_base - E2 idxs in baseline data
%frames_per_reward_range - a range on how many frames should elapse before a
% reward is expected.  Used to calibrate the target patterns.
%target_on_cov_bool - if false, calibrate target to baseline data directly.
% if true, use the covariance of data to calibrate target.  Potentially
% useful if the data itself does not exhibit strong co-modulation of E
% units.
%prefix_win - number of frames at start of baseline that calibration
% should NEGLECT.  (these starting frames are just ignored.)
%f0_win_bool - if true, estimate f0 with a window of activity.  if false, estimate f0 using the full baseline, 
%f0_win - number of frames to use for estimating f0.  (This code uses a
% rolling average, as used during the BMI)
%dff_win_bool - whether to smooth dff in a window
%dff_win - number of frames to use for smoothing dff
%save_dir - directory to save baseline results in in.

% %%
% %Debug parameters: 
% % (n_f_file, E_id, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
% 
% save_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline/results'
% mkdir(save_dir); 
% 
% target_on_cov_bool = 0; 
% 
% frame_rate = 10; 
% frames_per_reward_range = [1.5*60*frame_rate 1*60*frame_rate];
% prefix_win = 40; %number frames to remove from the start of imaging
% f0_win_bool = 1;
% f0_win = 2*60*frame_rate; %frames
% dff_win_bool = 1; 
% dff_win = 2; %frames
% 
% data_dir = '/Users/vivekathalye/Dropbox/Work/projects/holo_BMI/baseline'; 
% n_f_file = fullfile(data_dir, 'baseline_IntegrationRois_00001_1.csv'); 
% %Load baseline data:
% A = csvread(n_f_file, 1, 0);
% ts = A(:,1); 
% frame = A(:,2); 
% n = A(:, 3:end); 
% num_n = size(n,2); 
% [n_var, I] = sort(var(n, 0, 1), 'descend');
% n = n(:,I); 
% n(:,find(isnan(n_var))) = [];
% n_var(find(isnan(n_var))) = [];
% 
% direct = [1 2 4 5 6 17 27 39]; %Hand chosen for debugging
% f_raw = n(:,direct).'; %direct neuron fluorescence
% %Transposed to be num_neurons x num_samples, as expected for the input to
% %this file
% E_id = [2 2 2 2 1 1 1 1]'; 

%%
%Plots to show: 
%Plot colors: 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
%E1: blue-ish.  E2: orange-ish

plot_raw_bool = 1; 
plot_f0_bool = 1; 
plot_dff_bool = 1; 
plot_dff_z_smooth_bool = 1; 
plot_cov_bool = 1; 
%-

%%
%Select BMI ensemble temporal activity from baseline
%Create Ensemble Information
%Select BMI ensemble spatial components
%Decoder Information


%1) Select Temporal: BMI E data from baseline 

load(n_f_file); %num_samples x num_neurons_base
f_base = baseActivity.';
% f_base = f_base(1:1000, :); %debugging input data with nans... 
%Assume variable is called f_base

%Throw out prefix frames:
E1_raw = f_base((prefix_win+1):end,E1_base); 
E2_raw = f_base((prefix_win+1):end,E2_base); 
f_raw = [E1_raw E2_raw]; %first E1, then E2

%2) Ensemble information
num_E1 = length(E1_base); 
num_E2 = length(E2_base); 
num_neurons = num_E1 + num_E2;

E_id = [1*ones(num_E1, 1); 2*ones(num_E2, 1)]; 
E1_sel = E_id==1; 
E1_sel_idxs = find(E1_sel); 
E2_sel = E_id==2; 
E2_sel_idxs = find(E2_sel); 

%%
%3) Select Spatial Components for BMI E
% % Uncomment and test with a mosue:
load(Acomp_file, 'AComp');
E_base_sel = [E1_base, E2_base];
AComp_BMI = AComp(:, E_base_sel);
save(fullfile(save_dir, 'redcompBMI.mat'), 'AComp_BMI', 'E_base_sel', 'E_id');

%4) Decoder information
%Decoder information
E1_proj = zeros(num_neurons, 1); 
E1_proj(E1_sel) = 1;
E1_norm = sum(E1_sel); %can replace with vector norm.  

disp('E1 proj'); 
E1_proj = E1_proj/E1_norm;
E1_proj

E2_proj = zeros(num_neurons, 1); 
E2_proj(E2_sel) = 1; 
E2_norm = sum(E2_sel); 
disp('E2 proj')
E2_proj = E2_proj/E2_norm;
E2_proj

disp('decoder:')
decoder = E2_proj - E1_proj

%%
%--------------------------------------------------------------------------
%raw data plot: 
if(plot_raw_bool)
    [h, offset_vec] = plot_E_activity(f_raw, E_id, E_color);
    
%     h = figure;
%     hold on; 
%     offset = 0; 
%     for i=1:num_neurons
%         y_plot = f_raw(:,i);
%         y_plot = y_plot-min(y_plot); 
%         y_amp = max(y_plot); 
% 
%         if(i>1)
%             offset = offset+y_amp; 
%         end
% 
%         plot_color = E_color{E_id(i)};
%         plot(y_plot-offset, 'Color', plot_color); 
%     end
    xlabel('frame'); 
    ylabel('fluorescence'); 
    title('raw fluorescence in baseline'); 
    im_path = fullfile(save_dir, 'baseline_fraw.png'); 
    saveas(h, im_path); 
end

%%
%First process f0: 
if(f0_win_bool)
    %Implemented in the same fashion as in the BMI currently
    num_samples = size(f_raw,1);
    f0 = zeros(num_samples-f0_win+1, num_neurons); 
    f0(1,:) = mean(f_raw(1:f0_win, :), 1);
    for i = 2:length(f0)
        f0(i,:) = f0(i-1)*((f0_win-1)/f0_win) + f_raw((i+f0_win-1), :)/f0_win; 
    end
    %Truncate data based on the f0_win:
    f_postf0 = f_raw(f0_win:end, :); 
    f0_mean = repmat(mean(f_postf0, 1), size(f_postf0,1), 1);
else
    f_postf0 = f_raw; 
    f0_mean = repmat(mean(f_postf0, 1), size(f_postf0,1), 1);
    f0 = f0_mean; 
end

%%
%Compare f0win to f0mean:
if(plot_f0_bool)
    if(f0_win_bool)
        %Plot for one neuron: 
        n_i = 1; 
        h = figure; hold on;
        plot(f_postf0(:,n_i)); 
        plot(f0_mean(:,n_i), 'LineWidth', 5); 
        plot(f0(:,n_i), 'LineWidth', 5); 
        legend({'fraw', 'f0mean', 'f0win'})
%         h = figure; hold on; 
%         plot(f0(:,n_i)); 
%         plot(f0_mean(:,n_i)); 
%         legend({'f0win', 'f0mean'})
        xlabel('frame'); 
        ylabel('fluorescence'); 
        title('F0 for one neuron'); 
        saveas(h, fullfile(save_dir, 'f0.png')); 
    else
        %Plot for one neuron: 
        n_i = 1; 
        h = figure; hold on;
        plot(f_postf0(:,n_i)); 
        plot(f0_mean(:,n_i), 'LineWidth', 5); 
        legend({'fraw', 'f0mean'})
%         h = figure; hold on; 
%         plot(f0(:,n_i)); 
%         plot(f0_mean(:,n_i)); 
%         legend({'f0win', 'f0mean'})
        xlabel('frame'); 
        ylabel('fluorescence'); 
        title('F0 for one neuron');
        saveas(h, fullfile(save_dir, 'f0.png')); 
    end
end
%Note: a more sophisticated method would calculate f0 based on low-pass
%filtered calcium.  our f0 estimate is biased by ca transients.
%%
%Second, compute dff and dff_z:
dff = (f_postf0-f0)./f0;
%mean center the dff:
n_mean = mean(dff,1); %1 x num_neurons
mean_mat = repmat(n_mean, size(dff,1), 1);
dffc = dff-mean_mat;
%divide by std:
n_std = var(dffc, 0, 1).^(1/2); %1 x num_neurons
dff_z = dffc./repmat(n_std, [size(dff,1) 1]); 
if(plot_dff_bool)
    %plot dff
    h = plot_E_activity(dff, E_id, E_color);
    xlabel('frame'); 
    ylabel('dff'); 
    title('dff'); 
    im_path = fullfile(save_dir, 'dff.png'); 
    saveas(h, im_path);
    
    %plot dffz
    h = plot_E_activity(dff_z, E_id, E_color);
    xlabel('frame'); 
    ylabel('dff_z');    
    title('zscore dff'); 
    im_path = fullfile(save_dir, 'dffz.png'); 
    saveas(h, im_path); 
end

%%
%Third (final), smooth dff:
if(dff_win_bool)
    num_samples = max(length(dff_z)-(dff_win-1), 0); 
	n_analyze = zeros(num_samples, num_neurons); 
    smooth_filt = ones(dff_win,1)/dff_win; 
    for i=1:num_neurons
        n_analyze(:,i) = conv(dff_z(:,i), smooth_filt, 'valid'); 
    end
else
    n_analyze = dff_z; 
end

if(plot_dff_z_smooth_bool && dff_win_bool)
	h = figure; hold on;
    plot(dff_z(dff_win:end, 1)); 
    plot(n_analyze(:,1)); 
    legend({'dff z', 'dff z smooth'}); 
    xlabel('frame'); 
    ylabel('dff z'); 
    title('dff zscore vs smoothed dff zscore'); 
end

%%
analyze_cov = cov(n_analyze);
analyze_mean = mean(n_analyze); 
if(plot_cov_bool)
    h = figure;
    imagesc(analyze_cov); 
    colorbar
    axis square
    xlabel('roi')
    ylabel('roi')
    colormap;
    caxis([-0.2 0.3]); 
    title('neural cov'); 
    saveas(h, fullfile(save_dir, 'cov_mat_baseline.png'))

    [u,s,v] = svd(analyze_cov); 
    s_cumsum = cumsum(diag(s))/sum(diag(s)); 
    h = figure;
    plot(s_cumsum, '.-', 'MarkerSize', 20); 
    axis square
    xlabel('PC'); 
    ylabel('Frac Var Explained'); 
    title('DFF Z Smooth PCA Covariance');
    saveas(h, fullfile(save_dir, 'cov_pca_baseline.png'))
end

%%
%Cursor Cov:
cursor_cov = decoder'*analyze_cov*decoder; 
% cursor_cov = decoder'*test_cov*decoder; 
% cursor_cov

%%
%Inputs: 
%frames_per_reward_range
%cov_bool
reward_per_frame_range = 1./frames_per_reward_range;
E1_mean = mean(analyze_mean(E1_sel));
E1_std = sqrt((E1_sel/num_E1)'*analyze_cov*(E1_sel/num_E1));
E2_subord_mean = zeros(num_E2,1);
E2_subord_std = zeros(num_E2,1); 
E1_analyze = n_analyze(:,E1_sel); 
E2_analyze = n_analyze(:,E2_sel); 
for E2_i = 1:num_E2
    subord_sel = E2_sel;
    subord_sel(E2_sel_idxs(E2_i)) = 0; 
    E2_subord_mean(E2_i) = mean(analyze_mean(subord_sel));     
    var_i = subord_sel'*analyze_cov*subord_sel; 
    E2_subord_std(E2_i) = sqrt(var_i);     
end

E2_sum_analyze = sum(E2_analyze,2); 

%signals needed for target detection:
cursor_obs                      = n_analyze*decoder; 
E1_mean_analyze                 = mean(E1_analyze,2);
E1_mean_max                     = max(E1_mean_analyze); 
[E2_dom_samples, E2_dom_sel]    = max(E2_analyze, [], 2);
E2_subord_mean_analyze          = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);

%
% h = figure;
% hist(mean(E1_analyze,2)); 
% vline(E1_mean); 
% vline(E1_mean+E1_std); 

%%
%Iterate on T value, until perc correct value is achieved using truncated
%neural activity
T0 = max(n_analyze*decoder);
E2_coeff = 0.5; %multiples the std dev, for figuring out E2_subord_thresh
E2_subord_thresh = E2_subord_mean+E2_subord_std*E2_coeff;
E1_thresh = E1_mean + 2*E1_std; %E1_mean_max; %E1_mean;

T = T0; 
T_delta = 0.05;
E2_coeff_delta = 0.05; %0.05 
task_complete = 0;
T_vec = []; 
reward_per_frame_vec = []; 

max_iter = 10000;
iter = 0;

%If using data covariance:
rand_num_samples = 500000;
while(~task_complete)
    T_vec = [T_vec T];
    if(target_on_cov_bool)
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
        c2 = find(E1_mean_samples <= E1_thresh); 
        
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
        c2 = find(E1_mean_analyze <= E1_thresh);
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
        disp('target calibration complete!');
    elseif(reward_prob_per_frame > reward_per_frame_range(2))
        T = T+T_delta; 
    elseif(reward_prob_per_frame < reward_per_frame_range(1))
        %If we swept the full range of T, lower E2_coeff:
        if(T<=0)
            T=T0; 
            E2_coeff = E2_coeff - E2_coeff_delta;
            E2_subord_thresh = E2_subord_mean+E2_subord_std*E2_coeff; 
        end        
        T = T-T_delta; 
    end
    T
%     E2_coeff
%     E2_subord_thresh
    iter = iter+1;
    if(iter == max_iter)
        task_complete = 1;
        disp('Max Iter reached, check reward rate / baseline data'); 
    end
end

%%
h = figure;
plot(T_vec, '.-', 'MarkerSize', 7); 
xlabel('alg iteration'); 
ylabel('target'); 
title('Target Value over Calibration'); 
saveas(h, fullfile(save_dir, 'target_val_over_calibration.png')); 

% h = figure;
% hold on; 
% plot(reward_per_frame_vec, '.-', 'MarkerSize', 7); 
% hline(reward_per_frame_range(1)); 
% hline(reward_per_frame_range(2)); 
% xlabel('alg iteration'); 
% ylabel('reward per Frame'); 
% title('Reward Per Frame over Calibration'); 
% saveas(h, fullfile(save_dir, 'reward_per_frame_over_calibration.png')); 

%%
%Plot the hits on actual data: 
disp('T'); 
T
%data: 2.29
cursor_obs = n_analyze*decoder; 
c1 = find(cursor_obs >= T); 
disp('num E2-E1 >= T'); 
length(c1)

E1_mean_analyze = mean(E1_analyze,2)
c2 = find(E1_mean_analyze <= E1_thresh);
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
% vline(T); 
xlabel('Cursor'); 
ylabel('Number of Observations'); 
title('E2-E1 threshold plotted on E2-E1 distribution'); 
saveas(h, fullfile(save_dir, 'cursor_dist_T.png')); 

%%
%Plot the hit times: 
[h, offset_vec] = plot_E_activity(n_analyze, E_id, E_color);
xlabel('frame'); 
title(['Num Baseline Hits ' num2str(num_hits)]); 
offset = 5; 
%c1:
c1_offset = offset_vec(end)+offset;
plot(1:length(cursor_obs), cursor_obs-c1_offset);
% hline(T-c1_offset)

%c2:
c2_offset = offset_vec(end)+2*offset;
plot(1:length(E1_mean_analyze), E1_mean_analyze-c2_offset);
% hline(E1_thresh-c2_offset)

%c3:
c3_offset = offset_vec(end)+3*offset;
plot(E2_subord_mean_analyze-c3_offset); 
plot(E2_subord_thresh(E2_dom_sel)-c3_offset);

% for i=1:length(hit_times)
%     vline(hit_times(i)); 
% end

saveas(h, fullfile(save_dir, 'neural_hit_constraints.png')); 

%%
%Plot PSTH of neural activity locked to target hit: 
psth_win = [-30 30]*3; 
[psth_mean, psth_sem, psth_mat] = calc_psth(n_analyze, hit_times, psth_win);
h = figure; hold on;
offset = 0; 
for i=1:num_neurons
    y_plot = psth_mean(:,i); 
    y_plot = y_plot-min(y_plot);
    y_amp = max(y_plot); 
    offset = offset + y_amp; 
    y_sem = psth_sem(:,i)-min(y_plot); 
    
    plot(y_plot-offset, 'Color', E_color{(E_id(i))}); 
    errbar(1:length(y_plot), y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
end
% vline((psth_win(2)-psth_win(1))/2+1); 
xlabel('frame')
title('PSTH of Baseline Activity Locked to Target Hit'); 

saveas(h, fullfile(save_dir, 'PSTH_locked_to_hit_baseline.png')); 

% %%
% h = figure; hold on;
% for i =1:size(psth_mat,3)
%     plot(psth_mat(:,2,i)); 
% end
%
%%
%Save the results: 
%1) All the steps here
%2) Just the target parameters for running BMI

%1)All the steps here
clear h
date_str = datestr(datetime('now'), 'yyyymmddTHHMMSS'); 
save_path = fullfile(save_dir, ['target_calibration_ALL_' date_str '.mat']); 
save(save_path); 

%2)Just the target parameters for running BMI
save_path = fullfile(save_dir, 'BMI_target_info.mat'); 
%Change variable names for BMI code:
T1 = T; %Change to T1, as this is what BMI expects
save(save_path, 'n_mean', 'n_std', 'AComp_BMI', 'decoder', 'E_id', 'E1_sel_idxs', 'E2_sel_idxs', 'E1_base', 'E2_base', 'T1', 'E1_thresh', 'E2_subord_thresh', 'E2_coeff', 'E2_subord_mean', 'E2_subord_std'); 
