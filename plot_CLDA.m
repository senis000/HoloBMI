function [cal_update] = plot_CLDA(clda_mat, save_dir)

%%
%
%

%%
plot_bool = 1; 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 

%%
load(clda_mat); 
%data
%cal
%task_settings

%%
%Plot CLDA: 
valid_idxs      = ~isnan(data.E2_T); 
E2_T_valid      = data.E2_T(valid_idxs); 
E1_T_valid      = data.E1_T(valid_idxs); 
mid_T_valid     = data.mid_T(valid_idxs); 
h = figure;
plot(E2_T_valid); 
title('E2 T'); 

h = figure;
plot(E1_T_valid); 
title('E1 T'); 

h = figure;
plot(mid_T_valid); 
title('Mid T'); 
%%
%Get the correct number of hits with the final calibration parameters: 
cal_update = cal; 

E1_bool     = 1; 
E_cal       = cal.target.E1_hit_cal; 
cursor_obs  = data.cursor; 
mean_E2     = []; 
mean_E1     = []; 

b2base_coeff        = task_settings.b2base_coeff;
b2baseFrameThresh   = task_settings.back2BaseFrameThresh;

[E1_hit_cal, E1_hit_data] = compute_thresh_hits_b2base(E1_bool, E_cal, ...
    cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh);
cal_update.target.E1_hit_cal = E1_hit_cal; 

%%
E1_bool     = 0; 
E_cal       = cal.target.E2_hit_cal; 
cursor_obs  = data.cursor; 
mean_E2     = []; 
mean_E1     = []; 

b2base_coeff        = task_settings.b2base_coeff;
b2baseFrameThresh   = task_settings.back2BaseFrameThresh;

[E2_hit_cal, E2_hit_data] = compute_thresh_hits_b2base(E1_bool, E_cal, ...
    cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh);

cal_update.target.E2_hit_cal = E2_hit_cal; 
%%
%Summary results: 
disp('E2mE1 T:'); 
cal_update.target.E2_hit_cal.T

disp('E2 HIGH: num valid hits (WITH B2BASE):'); 
cal_update.target.E2_hit_cal.num_hits_b2base

disp('E2 HIGH: num baseline hits WITHOUT B2BASE:'); 
cal_update.target.E2_hit_cal.num_hits_no_b2base

disp('E1 HIGH: num valid hits (WITH B2BASE):'); 
cal_update.target.E1_hit_cal.num_hits_b2base

disp('E1 HIGH: num baseline hits WITHOUT B2BASE:'); 
cal_update.target.E1_hit_cal.num_hits_no_b2base

%%
%Plot cursor histogram: 

E2_T                = cal_update.target.E2_hit_cal.T;
E2_hits_b2base      = cal_update.target.E2_hit_cal.num_hits_b2base;
E2_hits_no_b2base   = cal_update.target.E2_hit_cal.num_hits_no_b2base;

E1_T                = cal_update.target.E1_hit_cal.T;
E1_hits_b2base      = cal_update.target.E1_hit_cal.num_hits_b2base;
E1_hits_no_b2base   = cal_update.target.E1_hit_cal.num_hits_no_b2base;

cursor_valid = data.cursor; 
cursor_valid(isnan(cursor_valid)) = []; 

h = figure;
hold on; 
hist(cursor_valid, 50); 
vline(E2_T); 
vline(E1_T); 

xlabel('Cursor'); 
ylabel('Number of Observations'); 

title(['E2: T: ' num2str(E2_T) ' hits b2base: ' num2str(E2_hits_b2base) ' hits no b2b: ' num2str(E2_hits_no_b2base) ...
    ' E1: T: ' num2str(E1_T) ' hits b2base: ' num2str(E1_hits_b2base) ' hits no b2b: ' num2str(E1_hits_no_b2base)]); 
saveas(h, fullfile(save_dir, 'cursor_dist_T.png'));


%%
%Plot auditory feedback
% plot_cursor = linspace(cal.fb.cursor_min, cal.fb.cursor_max, 1000); 
plot_cursor = linspace(min(data.cursor), max(data.cursor), 1000); 
plot_freq   = cursor2audio_freq_middle_match(plot_cursor, cal_update);
% plot_freq   = cursor2audio_freq(plot_cursor, cal);
h = figure;
hold on;
plot(plot_cursor, plot_freq); 
xlabel('Cursor E2-E1'); 
ylabel('Audiory Freq'); 
vline([cal_update.target.E1_hit_cal.T cal_update.target.E2_hit_cal.T]); 
% vline(); 
saveas(h, fullfile(save_dir, 'cursor2freq.png')); 

%%
% fb_obs = cursor2audio_freq(cursor_obs, cal);
cursor_obs_valid = cursor_obs; 
cursor_obs_valid(isnan(cursor_obs_valid)) = []; 


fb_obs = cursor2audio_freq_middle_match(cursor_obs_valid, cal);
num_fb_bins = 100; 
h = figure;
hist(fb_obs, num_fb_bins); 
xlabel('audio freq'); 
ylabel('baseline counts'); 
saveas(h, fullfile(save_dir, 'base_freq_hist.png')); 




%%
activity = data.bmiAct; 
activity(:, isnan(activity(1,:))) = []; 
activity = activity.'; %num_samples x num_neurons
num_neurons = size(activity, 2); 

%Throw out prefix frames, and re-order data:
E1_raw = activity((task_settings.prefix_win+1):end,cal.neurons.E1_sel_idxs);
E2_raw = activity((task_settings.prefix_win+1):end,cal.neurons.E2_sel_idxs); 
f_raw = [E1_raw E2_raw]; %first E1, then E2

%%
%First process f0: 
f0_win = task_settings.f0_win; 
if(task_settings.calibration.f0_win_bool)    
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

%%
E_id = cal_update.neurons.E_id;
[h, offset_vec] = plot_E_activity(1:length(f_postf0), f_postf0, E_id, E_color, 0);
xlabel('frame'); 
ylabel('fluorescence'); 
title('raw fluorescence in baseline'); 
im_path = fullfile(save_dir, 'baseline_fraw.png'); 
saveas(h, im_path); 

%%
%Second, smooth f:
num_samples = size(f_postf0,1);     
f_smooth = zeros(num_samples, num_neurons); 
smooth_filt = ones(task_settings.dff_win,1)/task_settings.dff_win;     
for i=1:num_neurons
    f_smooth(:,i) = conv(f_postf0(:,i), smooth_filt, 'same'); 
end

% if(plot_smooth_bool)
% 	h = figure; hold on;
%     plot(f_postf0(:,1)); 
%     plot(f_smooth(:,1)); 
%     legend({'f', 'f smooth'}); 
%     xlabel('frame'); 
%     ylabel('F'); 
%     title('F vs smoothed F'); 
% end

%%
%Third, compute dff and dff_z:
dff = (f_smooth-f0)./f0;
%mean center the dff:
n_mean = nanmean(dff,1); %1 x num_neurons
mean_mat = repmat(n_mean, size(dff,1), 1);
dffc = dff-mean_mat;
%divide by std:
n_std = nanstd(dffc, 0, 1); %var(dffc, 0, 1).^(1/2); %1 x num_neurons
dff_z = dffc./repmat(n_std, [size(dff,1) 1]); 

%range normalize: 
valid_idxs  = find(~isnan(dff(:,1))); 
dff_valid   = dff(valid_idxs, :);
n_max = prctile(dff_valid, 100, 1); 
n_min = prctile(dff_valid, 1, 1);
n_range = n_max - n_min; 
dff_min_c = dff-repmat(n_min, size(dff,1), 1);
dff_range_norm = dff_min_c./repmat(n_range, [size(dff,1) 1]); 
% dff_range_norm = dff./repmat(n_range, [size(dff,1) 1]); 

if(plot_bool)
    %plot dff
    h = plot_E_activity(1:length(dff), dff, E_id, E_color,0);
    xlabel('frame'); 
    ylabel('dff'); 
    title('dff'); 
    im_path = fullfile(save_dir, 'dff.png'); 
    saveas(h, im_path);
    
    %plot dffz
    if(task_settings.cursor_zscore_bool)
        h = plot_E_activity(1:length(dff_z), dff_z, E_id, E_color, 0);
        xlabel('frame'); 
        ylabel('dff_z');    
        title('zscore dff'); 
        im_path = fullfile(save_dir, 'dffz.png'); 
        saveas(h, im_path); 
    elseif(task_settings.range_norm_bool)
        h = plot_E_activity(1:length(dff_range_norm), dff_range_norm, E_id, E_color, 0);
        xlabel('frame'); 
        ylabel('dff_range_norm');    
        title('range norm dff'); 
        im_path = fullfile(save_dir, 'dffrangenorm.png'); 
        saveas(h, im_path);         
    end
end


%%
if task_settings.cursor_zscore_bool
    n_analyze = dff_z;
elseif task_settings.range_norm_bool
    n_analyze = dff_range_norm; 
else
    n_analyze = dff;
end
 
valid_idxs = find(~isnan(n_analyze(:,1)));
n_analyze = n_analyze(valid_idxs, :); 
analyze_cov = cov(n_analyze);
analyze_mean = nanmean(n_analyze); 
if(1)
    h = figure;
    imagesc(analyze_cov); 
    colorbar
    axis square
    xlabel('roi')
    ylabel('roi')
    colormap;
    caxis([-0.2 0.5]); 
    title('neural cov'); 
    saveas(h, fullfile(save_dir, 'cov_mat_baseline.png'))

    [u,s,v] = svd(analyze_cov); 
    s_cumsum = cumsum(diag(s))/sum(diag(s)); 
    h = figure;
    plot(s_cumsum, '.-', 'MarkerSize', 20); 
    axis square
    xlabel('PC'); 
    ylabel('Frac Var Explained'); 
    title('DFF Smooth PCA Covariance');
    saveas(h, fullfile(save_dir, 'cov_pca_baseline.png'))
end


%%
%Plot PSTH of neural activity locked to E2 target hit: 
valid_hit_idxs = E2_hit_data.T_idxs_b2base;

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

saveas(h, fullfile(save_dir, 'E2_PSTH_locked_to_hit_baseline.png')); 

%%
%Plot PSTH of neural activity locked to E2 target hit: 
valid_hit_idxs = E1_hit_data.T_idxs_b2base;

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

saveas(h, fullfile(save_dir, 'E1_PSTH_locked_to_hit_baseline.png')); 

