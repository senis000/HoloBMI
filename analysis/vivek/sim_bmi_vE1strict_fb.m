function [result_save_path] = sim_bmi_vE1strict_fb(n_data, task_settings, target_info, save_dir)

%%
% n_data = bmi_data.data.bmiAct; 
% sim_dir = fullfile(save_home, 'bmi_sim'); 
% save_dir = sim_dir

%%
%Input data / parameters: 
%n: num_neurons X num_samples

%%
plot_raw_bool = 1; 
plot_f0_bool = 1; 
plot_smooth_bool = 1; 
plot_dff_bool = 1; 
plot_cov_bool = 1; 

plotPath = fullfile(save_dir, 'plots'); 
mkdir(plotPath); 

%%
%1) remove nans: 
f = n_data; 


f(:,isnan(f(1,:))) = []; 
f = f.'; %num_samples x num_neurons

%2) Re-order data by E1, E2: 
if(size(f,2) > length(target_info.E_id))
    E1_temp = f(:,target_info.E1_base); 
    E2_temp = f(:,target_info.E2_base); 
    f = [E1_temp E2_temp];     
else
    E1_temp = f(:,target_info.E1_sel_idxs); 
    E2_temp = f(:,target_info.E2_sel_idxs); 
    f = [E1_temp E2_temp]; 
end
num_samples     = size(f,1); 
num_neurons     = size(f,2); 

%3) Throw out prefix frames:
prefix_win = task_settings.prefix_win;
f_raw = f((prefix_win+1):end,:); 

%4) Calculate f0
%f0_win_bool = 1; %can also load from task_settings: 
f0_win = task_settings.f0_win;
if task_settings.calibration.f0_win_bool
    %Calculate f0 as in BMI: 
    if task_settings.calibration.f0_init_slide
        num_samples = size(f_raw,1);
        f0 = zeros(num_samples,num_neurons); 
        for i=1:length(f0)
            if i==1
                f0(i,:) = f_raw(i,:);
            elseif i<f0_win
                f0(i,:) = (f0(i-1,:)*(i-1)+f_raw(i,:))/i;
            else
                f0(i,:) = (f0(i-1,:)*(f0_win-1)+f_raw(i,:))/f0_win;
            end
        end
        f_postf0 = f_raw;
        f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);
    else
        num_samples = size(f_raw,1);
        f0 = zeros(num_samples-f0_win+1, num_neurons); 
        f0(1,:) = mean(f_raw(1:f0_win, :), 1);
        for i = 2:length(f0)
            f0(i,:) = f0(i-1,:)*((f0_win-1)/f0_win) + f_raw((i+f0_win-1), :)/f0_win; 
        end
        %Truncate data based on the f0_win:
        f_postf0 = f_raw(f0_win:end, :); 
        f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);        
    end
else
    f_postf0 = f_raw; 
    f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);
    f0 = f0_mean; 
end

%5) Smooth
dff_win         = task_settings.dff_win; 
dff_win_bool    = 1; 
if(dff_win_bool)
    num_samples = size(f_postf0,1);     
	f_smooth = zeros(num_samples, num_neurons); 
    smooth_filt = ones(dff_win,1)/dff_win;     
    for i=1:num_neurons
        f_smooth(:,i) = conv(f_postf0(:,i), smooth_filt, 'same'); 
    end
else
    f_smooth = f_postf0; 
end

%6) dff and dff_z: 
dff = (f_smooth-f0)./f0;
%mean center the dff:
n_mean = nanmean(dff,1); %1 x num_neurons
mean_mat = repmat(n_mean, size(dff,1), 1);
dffc = dff-mean_mat;
%divide by std:
n_std = nanstd(dffc, 0, 1); %var(dffc, 0, 1).^(1/2); %1 x num_neurons
dff_z = dffc./repmat(n_std, [size(dff,1) 1]); 

cursor_zscore_bool = task_settings.cursor_zscore_bool;
if cursor_zscore_bool
    n_analyze = dff_z;
else
    n_analyze = dff;
end



%%
E1_sel              = target_info.E_id == 1; 
E1_sel_idxs         = target_info.E1_sel_idxs; 
num_E1              = length(E1_sel_idxs); 
E2_sel              = target_info.E_id == 2; 
E2_sel_idxs         = target_info.E2_sel_idxs; 
num_E2              = length(E2_sel_idxs); 

T                   = target_info.T1; 
E1_thresh           = target_info.E1_thresh; 
E2_subord_thresh    = target_info.E2_subord_thresh; 

%%
valid_idxs  = find(~isnan(n_analyze(:,1)));
n_analyze   = n_analyze(valid_idxs, :); 
analyze_cov = cov(n_analyze);
analyze_mean = nanmean(n_analyze); %takes mean along dim 1.  n_analyze is num_samples X num_neurons.



E1_mean                 = mean(analyze_mean(E1_sel));
E1_std                  = sqrt((E1_sel/num_E1)'*analyze_cov*(E1_sel/num_E1));
E2_subord_mean          = zeros(num_E2,1);
E2_subord_std           = zeros(num_E2,1); 
E1_analyze              = n_analyze(:,E1_sel); 
E2_analyze              = n_analyze(:,E2_sel); 
for E2_i = 1:num_E2
    subord_sel                      = E2_sel
    subord_sel(E2_sel_idxs(E2_i))   = 0; 
    E2_subord_mean(E2_i)            = mean(analyze_mean(subord_sel));     
    var_i                           = subord_sel'*analyze_cov*subord_sel; 
    E2_subord_std(E2_i)             = sqrt(var_i);     
end

E2_sum_analyze = sum(E2_analyze,2); 

%signals needed for target detection:
decoder                         = target_info.decoder;
cursor_obs                      = n_analyze*decoder; 
E1_mean_analyze                 = mean(E1_analyze,2);
E2_mean_analyze                 = mean(E2_analyze, 2); 
E1_mean_max                     = max(E1_mean_analyze); 
[E2_dom_samples, E2_dom_sel]    = max(E2_analyze, [], 2);
E2_subord_mean_analyze          = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);


%1) E2-E1 > alpha
c1 = find(cursor_obs >= T); 
%2) E1 < mu
c2 = find(E1_mean_analyze <= E1_thresh);
%3) E2_subord > mu (anded with previous constraint)
%For each idx, subtract the 
c3 = find(E2_subord_mean_analyze >= E2_subord_thresh(E2_dom_sel)); 
hit_idxs_no_b2base = intersect(intersect(c1, c2), c3);
%Remove hits that fall in a back2base

%----------------------------------------------------------------------
%Remove hits that fall in a back2base
b2base_thresh = 0.5*T;
hits_valid = ones(length(hit_idxs_no_b2base),1); 
back2BaseFramesThresh = task_settings.back2BaseFrameThresh; 
if length(hit_idxs_no_b2base) > 1
    for i = 2:length(hit_idxs_no_b2base)
        b2base_bool = sum(cursor_obs(hit_idxs_no_b2base(i-1):hit_idxs_no_b2base(i)) <= b2base_thresh) >= back2BaseFramesThresh; 
        hits_valid(i) = b2base_bool; 
    end
end
hit_idxs_b2base         = hit_idxs_no_b2base(find(hits_valid)); 
valid_hit_idxs          = hit_idxs_b2base;
reward_prob_per_frame   = sum(hits_valid)/length(n_analyze);    
%----------------------------------------------------------------------

%%
%Summary results:
disp('T'); 
T

% cursor_obs = n_analyze*decoder; 
% c1 = find(cursor_obs >= T); 
disp('num E2-E1 >= T'); 
num_c1 = length(c1)

% E1_mean_analyze = mean(E1_analyze,2)
% c2 = find(E1_mean_analyze <= E1_thresh);
disp('E1 >= b'); 
num_c2 = length(c2)

% E2_mean_analyze = mean(E2_analyze,2); 
% [E2_dom_samples, E2_dom_sel] = max(E2_analyze, [], 2);
% E2_subord_mean_analyze = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);
% %For each idx, subtract the 
% c3 = find(E2_subord_mean_analyze >= E2_subord_thresh(E2_dom_sel)); 
disp('E2 subord >= c'); 
num_c3 = length(c3)

num_cursor_hits = length(c1); 
disp('num cursor target hits (wo E1<thr, E2sub>thr :'); 
num_cursor_hits

disp('num baseline hits WITHOUT B2BASE:'); 
num_hits_no_b2base = length(hit_idxs_no_b2base)

disp('num valid hits (WITH B2BASE):'); 
num_valid_hits = length(valid_hit_idxs)

cursor_amp = (max(cursor_obs)-min(cursor_obs));
cursor_offset = cursor_amp/10; 
max_cursor = max(cursor_obs); 


% SAVE ALL
clear h
date_str = datestr(datetime('now'), 'yyyymmddTHHMMSS'); 
result_save_path = fullfile(save_dir, ['sim__ALL_' date_str '.mat']); 
save(result_save_path); 

%%
%PLOTS: 
%Plot colors: 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
%--------------------------------------------------------------------------
%raw data plot: 
E_id = target_info.E_id; 
if(plot_raw_bool)
    t_plot = 1:length(f_postf0);
    [h, offset_vec] = plot_E_activity(t_plot, f_postf0, E_id, E_color);
    
    xlabel('frame'); 
    ylabel('fluorescence'); 
    title('raw fluorescence in baseline'); 
    im_path = fullfile(plotPath, 'baseline_fraw.png'); 
    saveas(h, im_path); 
end

%%
%smooth f plot:

if(plot_smooth_bool && dff_win_bool)
	h = figure; hold on;
    plot(f_postf0(:,1)); 
    plot(f_smooth(:,1)); 
    legend({'f', 'f smooth'}); 
    xlabel('frame'); 
    ylabel('F'); 
    title('F vs smoothed F'); 
end

%%
%dff:
if(plot_dff_bool)
    %plot dff
    t_plot = 1:length(dff); 
    h = plot_E_activity(t_plot, dff, E_id, E_color);
    xlabel('frame'); 
    ylabel('dff'); 
    title('dff'); 
    im_path = fullfile(plotPath, 'dff.png'); 
    saveas(h, im_path);
    
    %plot dffz
    t_plot = 1:length(dff_z); 
    h = plot_E_activity(t_plot, dff_z, E_id, E_color);
    xlabel('frame'); 
    ylabel('dff_z');    
    title('zscore dff'); 
    im_path = fullfile(plotPath, 'dffz.png'); 
    saveas(h, im_path); 
end

%%
if(plot_cov_bool)
    h = figure;
    imagesc(analyze_cov); 
    colorbar
    axis square
    xlabel('roi')
    ylabel('roi')
    colormap;
    caxis([-0.2 0.5]); 
    title('neural cov'); 
    saveas(h, fullfile(plotPath, 'cov_mat_baseline.png'))

    [u,s,v] = svd(analyze_cov); 
    s_cumsum = cumsum(diag(s))/sum(diag(s)); 
    h = figure;
    plot(s_cumsum, '.-', 'MarkerSize', 20); 
    axis square
    xlabel('PC'); 
    ylabel('Frac Var Explained'); 
    title('DFF Smooth PCA Covariance');
    saveas(h, fullfile(plotPath, 'cov_pca_baseline.png'))
end

%%
h =figure; hold on;
scatter(c1, ones(length(c1),1)*max_cursor + cursor_offset, 15, 'r'); %plot(cursor_obs-cursor_offset, 'k'); 
scatter(c2, ones(length(c2),1)*max_cursor + 2*cursor_offset, 15, 'g'); %plot(cursor_obs-cursor_offset, 'k'); 
scatter(c3, ones(length(c3),1)*max_cursor + 3*cursor_offset, 15, 'b'); %plot(cursor_obs-cursor_offset, 'k'); 
plot(cursor_obs); 
hline(T); 
plot(E1_mean_analyze-cursor_amp); 
plot(E2_subord_mean_analyze-2*cursor_amp); 
xlabel('frame'); 
title(['hits with b2base: ' num2str(num_valid_hits)]); 
legend({'c1', 'c2 - E1 cond', 'c3 - E2 cond', 'cursor', 'E1 mean', 'E2 subord mean'}); 
vline(valid_hit_idxs); 
saveas(h, fullfile(plotPath, 'cursor_hit_ts.png')); 

%%
offset = 0; 
[h, offset_vec] = plot_cursor_E1_E2_activity(cursor_obs, E1_mean_analyze, E2_mean_analyze, n_analyze, E_id, E_color, offset)
hold on; hline(T); 
saveas(h, fullfile(plotPath, 'cursor_E1_E2_ts.png')); 

%%
h = figure;
hold on; 
hist(cursor_obs, 50); 
vline(T); 
xlabel('Cursor'); 
ylabel('Number of Observations'); 
title(['E2-E1 thr on E2-E1 hist, num valid hits: ' num2str(num_valid_hits) ...
    ' num hits no b2base: ' num2str(num_hits_no_b2base) ...
    ' num cursor hits: ' num2str(num_cursor_hits)]); 
saveas(h, fullfile(plotPath, 'cursor_dist_T.png')); 

%%
%Plot PSTH of neural activity locked to target hit: 
psth_win = [-30 30]*3; 
[psth_mean, psth_sem, psth_mat] = calc_psth(n_analyze, valid_hit_idxs, psth_win);
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

saveas(h, fullfile(plotPath, 'PSTH_locked_to_hit_baseline.png')); 