function [cal, BMI_roi_path] = baseline2two_target_linear_fb(n_f_file, roi_data_file, task_settings, ...
    E1_base, E2_base, save_dir)


%%
%Initialize parameters to calibrate: 
%paths
%--------------------------------------------------------------------------
cal.paths.BMI_roi_path = ...
    fullfile(save_dir, 'BMI_roi.mat');
plotPath = fullfile(save_dir, 'plots');
cal.paths.plot_path = ...
    plotPath;
date_str            = datestr(datetime('now'), 'yyyymmddTHHMMSS'); 
cal_all_path        = fullfile(save_dir, ['BMI_cal_ALL_' date_str '.mat']); 
cal.paths.cal_all   = cal_all_path;
cal_path            = fullfile(save_dir, ['BMI_cal_' date_str '.mat']);
cal.paths.cal       = cal_path;

%reward
%--------------------------------------------------------------------------
cal.reward.sec_per_reward_range = ...
    task_settings.calibration.sec_per_reward_range;
cal.reward.num_base_samples = ...
    []; %assigned
cal.reward.baseline_frameRate = ...
    []; %assigned
cal.reward.frames_per_reward_range = ...
    []; %assigned

%BMI neurons
%--------------------------------------------------------------------------
[E_base_sel, num_E1, num_E2, num_neurons, ...
    E_id, E1_sel, E2_sel, E1_sel_idxs, E2_sel_idxs] = ...
    ensemble_info(E1_base, E2_base);
cal.neurons.strcMask = ...
    []; %assigned
cal.neurons.E1_base = ...
    E1_base;
cal.neurons.E2_base = ...
    E2_base;
cal.neurons.E_base_sel = ...
    E_base_sel;
cal.neurons.num_E1 = ...
    num_E1; 
cal.neurons.num_E2 = ...
    num_E2; 
cal.neurons.num_neurons = ...
    num_neurons; 
cal.neurons.E_id = ...
    E_id; 
cal.neurons.E1_sel = ...
    E1_sel; 
cal.neurons.E2_sel = ...
    E2_sel; 
cal.neurons.E1_sel_idxs = ...
	E1_sel_idxs; 
cal.neurons.E2_sel_idxs = ...
    E2_sel_idxs;  

%Decoder
%--------------------------------------------------------------------------
[decoder, E1_proj, E2_proj, E1_norm, E2_norm] = ...
    def_decoder(num_neurons, E1_sel, E2_sel);
cal.decoder = decoder; 
cal.cursor_zscore_bool = task_settings.cursor_zscore_bool; 

%Target (assigned later)
%--------------------------------------------------------------------------
cal.target.T = ...
    []; 
cal.target.E1_thresh = ...
    []; 
cal.target.E2_subord_thresh = ...
    [];
cal.target.n_mean = ...
    [];
cal.target.n_std = ...
    [];
cal.target.E1_mean = ...
    [];
cal.target.E1_std = ...
    [];
cal.target.E1_coeff = ...
    [];
cal.target.E2_subord_mean = ...
    [];
cal.target.E2_subord_std = ...
    [];
cal.target.E2_coeff = ...
    [];

% %Result (assigned later)
% %--------------------------------------------------------------------------
% cal.result.num_c1 = ...
%     []; 
% cal.result.num_c2 = ...
%     []; 
% cal.result.num_c3 = ...
%     []; 
% cal.result.num_hits_no_b2base = ...
%     []; 
% cal.result.num_valid_hits = ...
%     []; 

% %Auditory feedback 
% %--------------------------------------------------------------------------
% cal.fb.target_low_freq = ...
%     task_settings.fb.target_low_freq; 
% cal.fb.freq_min = ...
%     task_settings.fb.freq_min; 
% cal.fb.freq_min = ...
%     task_settings.fb.freq_min; 
% %Fields updated in calibration: 
% cal.fb.cursor_min         = ...
%     []; %for fb, cursor is ceil to this value
% cal.fb.cursor_max         = ...
%     []; %for fb, cursor is floor to this value
% cal.fb.cursor_range       = ...
%     []; 
% % freq = a*exp(b*(cursor_trunc-cursor_min))
% cal.fb.a                  = ...
%     []; 
% cal.fb.b                  = ...
%     []; 

%Variables formerly saved for the BMI: 
% neuron:
% 'AComp_BMI' , now: strcMask
% 'E_id', 'E1_sel_idxs', 'E2_sel_idxs', 'E1_base', 'E2_base'
%target:
%  'n_mean', 'n_std'
%'T1', 'E1_thresh', 'E1_coeff', 'E1_std', 'E2_subord_thresh', 'E2_coeff', 'E2_subord_mean', 'E2_subord_std' 
%decoder:
% 'decoder'

%%
%Plots to show: 
%Plot colors: 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
%E1: blue-ish.  E2: orange-ish

plot_bool           = 0; 
plot_raw_bool       = 1 && plot_bool; 
plot_f0_bool        = 1 && plot_bool; 
plot_smooth_bool    = 1 && plot_bool; 
plot_dff_bool       = 1 && plot_bool; 
plot_cov_bool       = 1 && plot_bool; 
%-

plotPath = fullfile(save_dir, 'plots'); 
mkdir(plotPath); 
%%
%Select BMI ensemble temporal activity from baseline
%Create Ensemble Information
%Select BMI ensemble spatial components
%Decoder Information


%1) Select Temporal: BMI E data from baseline 

load(n_f_file); 
f_base = baseActivity;  %num_neurons X num_samples
f_base(:,isnan(f_base(1,:))) = []; %remove nan data
f_base = f_base.'; %num_samples X num_neurons

%%
E1_temp = f_base(:,E1_base); 
E2_temp = f_base(:,E2_base); 
f = [E1_temp E2_temp]; 

%Throw out prefix frames:
E1_raw = f_base((task_settings.prefix_win+1):end,E1_base); 
E2_raw = f_base((task_settings.prefix_win+1):end,E2_base); 
f_raw = [E1_raw E2_raw]; %first E1, then E2

%%
%3) Select Spatial Components for BMI E
% % Uncomment and test with a mosue:

if task_settings.onacid_bool
    load(task_settings.Acomp_file, 'AComp');
    AComp_BMI = AComp(:, E_base_sel);
    strcMask = obtainStrcMask(AComp_BMI, px, py);
    save(fullfile(save_dir, 'redcompBMI.mat'), 'strcMask', 'E_base_sel', 'E_id');
else
    AComp_BMI = [];
    load(roi_data_file, 'roi_data'); 
    roi_mask = roi_data.roi_mask;
    EnsembleMask = zeros(size(roi_mask));
    for indn = 1:length(E1_base)
        auxmask = roi_mask;
        auxmask(auxmask~=E1_base(indn)) = 0;
        auxmask(auxmask~=0) = indn;
        EnsembleMask = auxmask + EnsembleMask;
    end
    for indn = 1:length(E2_base)
        auxmask = roi_mask;
        auxmask(auxmask~=E2_base(indn)) = 0;
        auxmask(auxmask~=0) = indn + length(E1_base);
        EnsembleMask = auxmask + EnsembleMask;
    end
    strcMask = obtainStrcMaskfromMask(EnsembleMask);
    BMI_roi_path = fullfile(save_dir, 'BMI_roi.mat');
    save(BMI_roi_path, 'strcMask', 'E_base_sel', 'E_id'); 
    %ASSIGN:
    cal.neurons.strcMask = ...
        strcMask;
    %In old version: 
%     save(fullfile(save_dir, 'redcompBMI.mat'), 'strcMask', 'E_base_sel', 'E_id'); 
end

%%
%4) Decoder information
%Decoder information
[decoder, E1_proj, E2_proj, E1_norm, E2_norm] = ...
    def_decoder(num_neurons, E1_sel, E2_sel); 

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
%Convert sec_per_reward_range to frames_per_reward_range: 
%--------------------------------------------------------------------------
cal.reward.num_base_samples = ...
    size(f_postf0,1); 
cal.reward_baseline_len_usable = ...
    size(f_postf0,1)/size(f_base,1)*task_settings.calibration.baseline_len;

%baseline_frameRate
cal.reward.baseline_frameRate = ...
    cal.reward.num_base_samples/cal.reward_baseline_len_usable;
%frames_per_reward_range
cal.reward.frames_per_reward_range = ...
    task_settings.calibration.sec_per_reward_range*cal.reward.baseline_frameRate; 

cal.reward.reward_per_frame_range = ...
    1./cal.reward.frames_per_reward_range;
%%
%--------------------------------------------------------------------------
%raw data plot: 
if(plot_raw_bool)
    [h, offset_vec] = plot_E_activity(1:length(f_postf0), f_postf0, E_id, E_color, 0);
%      plot_E_activity(t,n, E_id, E_color, offset)
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
    im_path = fullfile(plotPath, 'baseline_fraw.png'); 
    saveas(h, im_path); 
end

%%
%Compare f0win to f0mean:
if(plot_f0_bool)
    if(task_settings.calibration.f0_win_bool)
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
        saveas(h, fullfile(plotPath, 'f0.png')); 
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
        saveas(h, fullfile(plotPath, 'f0.png')); 
    end
end
%Note: a more sophisticated method would calculate f0 based on low-pass
%filtered calcium.  our f0 estimate is biased by ca transients.

%%
%Second, smooth f:
num_samples = size(f_postf0,1);     
f_smooth = zeros(num_samples, num_neurons); 
smooth_filt = ones(task_settings.dff_win,1)/task_settings.dff_win;     
for i=1:num_neurons
    f_smooth(:,i) = conv(f_postf0(:,i), smooth_filt, 'same'); 
end

if(plot_smooth_bool)
	h = figure; hold on;
    plot(f_postf0(:,1)); 
    plot(f_smooth(:,1)); 
    legend({'f', 'f smooth'}); 
    xlabel('frame'); 
    ylabel('F'); 
    title('F vs smoothed F'); 
end

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
if(plot_dff_bool)
    %plot dff
    h = plot_E_activity(1:length(dff), dff, E_id, E_color,0);
    xlabel('frame'); 
    ylabel('dff'); 
    title('dff'); 
    im_path = fullfile(plotPath, 'dff.png'); 
    saveas(h, im_path);
    
    %plot dffz
    h = plot_E_activity(1:length(dff_z), dff_z, E_id, E_color, 0);
    xlabel('frame'); 
    ylabel('dff_z');    
    title('zscore dff'); 
    im_path = fullfile(plotPath, 'dffz.png'); 
    saveas(h, im_path); 
end


%%
if task_settings.cursor_zscore_bool
    n_analyze = dff_z;
else
    n_analyze = dff;
end
 
valid_idxs = find(~isnan(n_analyze(:,1)));
n_analyze = n_analyze(valid_idxs, :); 
analyze_cov = cov(n_analyze);
analyze_mean = nanmean(n_analyze); 
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
%Cursor Cov:

cursor_cov = decoder'*analyze_cov*decoder; 
% cursor_cov = decoder'*test_cov*decoder; 
% cursor_cov

%%
%Inputs: 
%frames_per_reward_range
%cov_bool
reward_per_frame_range = 1./cal.reward.frames_per_reward_range
E1_mean = mean(analyze_mean(E1_sel));
E1_std = sqrt((E1_sel/num_E1)'*analyze_cov*(E1_sel/num_E1));

E2_mean = mean(analyze_mean(E2_sel));
E2_std = sqrt((E2_sel/num_E2)'*analyze_cov*(E2_sel/num_E2));

E1_analyze = n_analyze(:,E1_sel); 
E2_analyze = n_analyze(:,E2_sel);

E1_mean_analyze = mean(E1_analyze, 2); 
E2_mean_analyze = mean(E2_analyze, 2); 

%signals needed for target detection:
cursor_obs                      = n_analyze*decoder;
E2mE1_obs                       = cursor_obs; 
E1mE2_obs                       = -cursor_obs; 
%
% h = figure;
% hist(mean(E1_analyze,2)); 
% vline(E1_mean); 
% vline(E1_mean+E1_std); 

b2base_coeff = task_settings.b2base_coeff; 
b2baseFrameThresh = task_settings.back2BaseFrameThresh; 

%T:
T_high  = task_settings.fb.max_prctile; 
T_low   = task_settings.fb.min_prctile; 

T0      = task_settings.fb.max_prctile; 
T       = T0;
T_delta_min = 1e-7; 

%E2:
E2_hit_cal.E2_coeff                 = 0.0; %multiples the std dev
E2_hit_cal.E2_thresh                = E2_mean + E2_hit_cal.E2_coeff*E2_std;
E2_hit_cal.E1_coeff                 = 0.0; %multiples the std dev
E2_hit_cal.E1_thresh                = E1_mean - E2_hit_cal.E1_coeff*E1_std; 

%E1:
E1_hit_cal.E1_coeff                 = 0.0; 
E1_hit_cal.E1_thresh                = E1_mean + E1_hit_cal.E1_coeff*E1_std;
E1_hit_cal.E2_coeff                 = 0.0; 
E1_hit_cal.E2_thresh                = E2_mean + E1_hit_cal.E2_coeff*E2_std;

%Step size: 
task_complete               = 0;
E2_complete                 = 0;
E1_complete                 = 0; 
abort_calibration           = 0; 
T_vec                       = []; 
E1_reward_per_frame_vec     = []; 
E2_reward_per_frame_vec     = [];

init_best_cal               = 1; 
best_cal.n_mean             = n_mean;
best_cal.n_std              = n_std;
best_cal.E2_mean            = E2_mean; 
best_cal.E1_mean            = E1_mean; 
best_cal.E2_std             = E2_std; 
best_cal.E1_std             = E1_std; 

best_cal.E1_hit_cal             = E1_hit_cal; 
best_cal.E2_hit_cal             = E2_hit_cal; 

best_cal.T_prctile          = []; 
best_cal.score              = 1000; 
best_cal_data.E1_hit_data         = []; 
best_cal_data.E2_hit_data         = []; 

max_iter    = 1000;
iter        = 0;

%1) Find an E2 calibration which works
while(~task_complete)
    iter
    if(iter > max_iter)
        break;
    end
    
    T_vec = [T_vec T]; 
    T_prctile = T;
    cursor_analyze = cursor_obs;
    
    E1_bool = 0; 
    [E2_hit_cal, E2_hit_data] = ...
        compute_hits_b2base(E1_bool, E2_hit_cal, ...
        T_prctile, cursor_analyze, ...
        E2_mean_analyze, E1_mean_analyze, ...
        b2base_coeff, b2baseFrameThresh);
    
    E1_bool = 1; 
    E1_prctile = 100-T_prctile;
    [E1_hit_cal, E1_hit_data] = ...
        compute_hits_b2base(E1_bool, E1_hit_cal, ...
        E1_prctile, cursor_analyze, ...
        E2_mean_analyze, E1_mean_analyze, ...
        b2base_coeff, b2baseFrameThresh);      
    
    %Return Values:
    %T_value
    %T_idxs_no_b2base
    %b2base_thresh
    %T_idxs_b2base
    %num_hits_b2base
    %num_hits_no_b2base
    %reward_prob_per_frame

    E2_reward_prob = E2_hit_data.reward_prob_per_frame;
    E1_reward_prob = E1_hit_data.reward_prob_per_frame;
    [score, E2_score, E1_score] = ...
        score_calibration_2target(E1_reward_prob, E2_reward_prob, reward_per_frame_range);

    update_best_cal = ((E2_score == 0) && score <= best_cal.score); 
    if(update_best_cal)
        init_best_cal = 0; 
        T
        score
        disp('update best cal'); 
        best_cal.T_prctile      = T; 
        best_cal.score          = score;
        best_cal.E1_hit_cal     = E1_hit_cal; 
        best_cal.E2_hit_cal     = E2_hit_cal; 
        best_cal_data.E1_hit_data     = E1_hit_data; 
        best_cal_data.E2_hit_data     = E2_hit_data;         
    end

    if(E2_score == 0 && E1_score == 0)
        disp('task complete!'); 
        task_complete = 1; 
    else
%         reward_prob_per_frame = E2_reward_prob; 
        if(E2_score ~=0)
            reward_prob_per_frame = E2_reward_prob; 
        else
            disp('optimizing E1'); 
            reward_prob_per_frame = E1_reward_prob; 
        end

        if(reward_prob_per_frame > reward_per_frame_range(2))
            %Task too easy, make T harder:
            T_low = T
            T = (T_high + T_low)/2;
        elseif(reward_prob_per_frame < reward_per_frame_range(1))
            %Task too hard, make T easier:
            T_high = T
            T = (T_high +T_low)/2;
        end
    end
    if(abs(T-T_vec(end)) <= T_delta_min)
        T
        T_vec(end)
        T-T_vec(end)
        break;
    end
    iter = iter+1; 
end

%%
%Target parameters: 
cal.target = best_cal;

%%
%Summary results: 
disp('E2mE1 T:'); 
best_cal.E2_hit_cal.T

disp('E2 HIGH: num valid hits (WITH B2BASE):'); 
best_cal.E2_hit_cal.num_hits_b2base

disp('E2 HIGH: num baseline hits WITHOUT B2BASE:'); 
best_cal.E2_hit_cal.num_hits_no_b2base

disp('E1 HIGH: num valid hits (WITH B2BASE):'); 
best_cal.E1_hit_cal.num_hits_b2base

disp('E1 HIGH: num baseline hits WITHOUT B2BASE:'); 
best_cal.E1_hit_cal.num_hits_no_b2base

%%
% Calculate parameters for auditory feedback
[cal]   = cursor2audio_fb(cal, cursor_obs, task_settings);

%%
%PLOTS
%TODO: decomposition into plotting functions
%--------------------------------------------------------------------------
%%
%Plot auditory feedback
% plot_cursor = linspace(cal.fb.cursor_min, cal.fb.cursor_max, 1000); 
plot_cursor = linspace(min(cursor_obs), max(cursor_obs), 1000); 
plot_freq   = cursor2audio_freq(plot_cursor, cal);
h = figure;
hold on;
plot(plot_cursor, plot_freq); 
xlabel('Cursor E2-E1'); 
ylabel('Audiory Freq'); 
vline([best_cal.E1_hit_cal.T best_cal.E2_hit_cal.T]); 
% vline(); 
saveas(h, fullfile(plotPath, 'cursor2freq.png')); 

%%
fb_obs = cursor2audio_freq(cursor_obs, cal);
num_fb_bins = 100; 
h = figure;
hist(fb_obs, num_fb_bins); 
xlabel('audio freq'); 
ylabel('baseline counts'); 
saveas(h, fullfile(plotPath, 'base_freq_hist.png')); 

%%
plot_cal_hits(best_cal_data.E2_hit_data,...
    E1_mean_analyze, E2_mean_analyze, ...
    cursor_obs, 'E2 hits', fullfile(plotPath, 'E2_hit_ts.png')); 
%%
plot_cal_hits(best_cal_data.E1_hit_data,...
    E1_mean_analyze, E2_mean_analyze, ...
    cursor_obs, 'E1 hits', fullfile(plotPath, 'E1_hit_ts.png')); 
%%
offset = 0; 
[h, offset_vec] = plot_cursor_E1_E2_activity(cursor_obs, E1_mean_analyze, E2_mean_analyze, n_analyze, E_id, E_color, offset)
hold on; hline(cal.target.E2_hit_cal.T); 
saveas(h, fullfile(plotPath, 'cursor_E1_E2_ts.png')); 

%%
% hit_cal
% hit_cal.T                   = T_value;
% hit_cal.num_E2_valid        = num_E2_valid; 
% hit_cal.num_E1_valid        = num_E1_valid; 
% hit_cal.num_hits_b2base     = num_hits_b2base; 
% hit_cal.num_hits_no_b2base  = num_hits_no_b2base; 
% hit_cal.reward_prob_per_frame   = reward_prob_per_frame; 
%
% hit_data
% hit_data                    = hit_cal;
% hit_data.T_hits             = T_hits; 
% hit_data.E2_valid           = E2_valid; 
% hit_data.E1_valid           = E1_valid; 
% hit_data.T_idxs_no_b2base   = T_idxs_no_b2base;
% hit_data.T_idxs_b2base      = T_idxs_b2base;

E2_T                = cal.target.E2_hit_cal.T;
E2_hits_b2base      = cal.target.E2_hit_cal.num_hits_b2base;
E2_hits_no_b2base   = cal.target.E2_hit_cal.num_hits_no_b2base;

E1_T                = cal.target.E1_hit_cal.T;
E1_hits_b2base      = cal.target.E1_hit_cal.num_hits_b2base;
E1_hits_no_b2base   = cal.target.E1_hit_cal.num_hits_no_b2base;

h = figure;
hold on; 
hist(cursor_obs, 50); 
vline(E2_T); 
vline(E1_T); 

xlabel('Cursor'); 
ylabel('Number of Observations'); 

title(['E2: T: ' num2str(E2_T) ' hits b2base: ' num2str(E2_hits_b2base) ' hits no b2b: ' num2str(E2_hits_no_b2base) ...
    ' E1: T: ' num2str(E1_T) ' hits b2base: ' num2str(E1_hits_b2base) ' hits no b2b: ' num2str(E1_hits_no_b2base)]); 
saveas(h, fullfile(plotPath, 'cursor_dist_T.png'));

% %%
% %Plot PSTH of neural activity locked to E2 target hit: 
% valid_hit_idxs = best_cal_data.E2_hit_data.T_idxs_b2base;
% 
% psth_win = [-30 30]*3; 
% [psth_mean, psth_sem, psth_mat] = calc_psth(n_analyze, valid_hit_idxs, psth_win);
% h = figure; hold on;
% offset = 0; 
% for i=1:num_neurons
%     y_plot = psth_mean(:,i); 
%     y_plot = y_plot-min(y_plot);
%     y_amp = max(y_plot); 
%     offset = offset + y_amp; 
%     y_sem = psth_sem(:,i)-min(y_plot); 
%     
%     plot(y_plot-offset, 'Color', E_color{(E_id(i))}); 
%     errbar(1:length(y_plot), y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
% end
% % vline((psth_win(2)-psth_win(1))/2+1); 
% xlabel('frame')
% title('PSTH of Baseline Activity Locked to Target Hit'); 
% 
% saveas(h, fullfile(plotPath, 'E2_PSTH_locked_to_hit_baseline.png')); 
% 
% %%
% %Plot PSTH of neural activity locked to E2 target hit: 
% valid_hit_idxs = best_cal_data.E1_hit_data.T_idxs_b2base;
% 
% psth_win = [-30 30]*3; 
% [psth_mean, psth_sem, psth_mat] = calc_psth(n_analyze, valid_hit_idxs, psth_win);
% h = figure; hold on;
% offset = 0; 
% for i=1:num_neurons
%     y_plot = psth_mean(:,i); 
%     y_plot = y_plot-min(y_plot);
%     y_amp = max(y_plot); 
%     offset = offset + y_amp; 
%     y_sem = psth_sem(:,i)-min(y_plot); 
%     
%     plot(y_plot-offset, 'Color', E_color{(E_id(i))}); 
%     errbar(1:length(y_plot), y_plot-offset,y_sem, 'Color', E_color{(E_id(i))}); 
% end
% % vline((psth_win(2)-psth_win(1))/2+1); 
% xlabel('frame')
% title('PSTH of Baseline Activity Locked to Target Hit'); 
% 
% saveas(h, fullfile(plotPath, 'E1_PSTH_locked_to_hit_baseline.png')); 

%%
%Save the results: 
%1) All the steps here
%2) Just the target parameters for running BMI

%1)All the steps here
clear h
save(cal.paths.cal_all); 

%2)Just the calibration parameters for running BMI
save(cal.paths.cal, 'cal'); 

end

function [E_base_sel, num_E1, num_E2, num_neurons, ...
    E_id, E1_sel, E2_sel, E1_sel_idxs, E2_sel_idxs] = ...
    ensemble_info(E1_base, E2_base)

E_base_sel = [E1_base, E2_base]; 
num_E1 = length(E1_base); 
num_E2 = length(E2_base); 
num_neurons = ...
    num_E1 + num_E2;  
E_id = ...
    [1*ones(num_E1, 1); 2*ones(num_E2, 1)];
E1_sel = ...
    E_id == 1;
E2_sel = ...
    E_id == 2;
E1_sel_idxs = ...
	find(E1_sel); 
E2_sel_idxs = ...
    find(E2_sel); 
end

function [decoder, E1_proj, E2_proj, E1_norm, E2_norm] = ...
    def_decoder(num_neurons, E1_sel, E2_sel)

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
end

function [cal] = cursor2audio_fb(cal, cursor_obs, task_settings)
%INPUT: 
%cursor_obs
%cal
%FIELDS: 
% -- target.T
% -- fb.freq_min
% -- fb.fre_max

%Copy from task_settings:
cal.fb.fb_bool = ...
    task_settings.fb.fb_bool; 
cal.fb.target_low_freq = ...
    task_settings.fb.target_low_freq; 
cal.fb.freq_min = ...
    task_settings.fb.freq_min; 
cal.fb.freq_max = ...
    task_settings.fb.freq_max; 

%Calculate:
% cal.fb.min_perctile     = ...
%     10; 
cal.fb.cursor_min         = ...
    cal.target.E1_hit_cal.T; 
%negate because E1 was calibrated to E1mE2

cal.fb.cursor_max         = ...
    cal.target.E2_hit_cal.T; %for fb, cursor is floor to this value
cal.fb.cursor_range       = ...
    cal.fb.cursor_max - cal.fb.cursor_min; 
% freq = a*exp(b*(cursor_trunc-cursor_min))
cal.fb.a                  = ...
    cal.fb.freq_min; 
cal.fb.b                  = ...
    (log(cal.fb.freq_max) - log(cal.fb.a))/cal.fb.cursor_range; 

end


% compute_hits_b2base(T_prctile, cursor_analyze, ...
%         E1_mean_analyze, E1_thresh_high, ...
%         E2_mean_analyze, E2_thresh_low, ...
%         b2base_coeff, b2baseFrameThresh);


function [E_hit_cal, E_hit_data] = compute_hits_b2base(E1_bool, E_cal, ...
    T_prctile, cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh)
%Return Values:
% hit_cal
% hit_cal.T                   = T_value;
% hit_cal.num_E2_valid        = num_E2_valid; 
% hit_cal.num_E1_valid        = num_E1_valid; 
% hit_cal.num_hits_b2base     = num_hits_b2base; 
% hit_cal.num_hits_no_b2base  = num_hits_no_b2base; 
% hit_cal.reward_prob_per_frame   = reward_prob_per_frame; 
%
% hit_data
% hit_data                    = hit_cal;
% hit_data.T_hits             = T_hits; 
% hit_data.E2_valid           = E2_valid; 
% hit_data.E1_valid           = E1_valid; 
% hit_data.T_idxs_no_b2base   = T_idxs_no_b2base;
% hit_data.T_idxs_b2base      = T_idxs_b2base;


T_value                         = prctile(cursor_obs, T_prctile); 
b2base_thresh                   = b2base_coeff*T_value;

if(E1_bool)
%     T_neg                           = prctile(-cursor_obs, T_prctile); 
%     T_value                         = -T_neg; 
%     b2base_thresh                   = b2base_coeff*T_value;        
    T_hits                          = find(cursor_obs <= T_value); 
    E2_valid                        = find(mean_E2 <= E_cal.E2_thresh); 
    E1_valid                        = find(mean_E1 >= E_cal.E1_thresh);    
else
    T_hits                          = find(cursor_obs >= T_value); 
    E2_valid                        = find(mean_E2 >= E_cal.E2_thresh); 
    E1_valid                        = find(mean_E1 <= E_cal.E1_thresh);
end

%Intersection of these are the target hits without b2base constraint:
T_idxs_no_b2base        = intersect(intersect(E1_valid, E2_valid), T_hits); 
 
hits_valid              = ones(length(T_idxs_no_b2base),1); 
if length(T_idxs_no_b2base) > 1
    for i = 2:length(T_idxs_no_b2base)
        if(E1_bool)
            b2base_bool = sum(cursor_obs(T_idxs_no_b2base(i-1):T_idxs_no_b2base(i)) >= b2base_thresh) >= b2baseFrameThresh; 
        else
            b2base_bool = sum(cursor_obs(T_idxs_no_b2base(i-1):T_idxs_no_b2base(i)) <= b2base_thresh) >= b2baseFrameThresh; 
        end
        hits_valid(i) = b2base_bool; 
    end
end
T_idxs_b2base               = T_idxs_no_b2base(find(hits_valid));    

num_E2_valid                = length(E2_valid); 
num_E1_valid                = length(E1_valid); 
num_hits_b2base             = length(T_idxs_b2base); 
num_hits_no_b2base          = length(T_idxs_no_b2base); 
reward_prob_per_frame       = ...
    num_hits_b2base/length(cursor_obs);

%ASSIGN:
%hit_cal
E_hit_cal                       = E_cal; 
E_hit_cal.T                   = T_value;
E_hit_cal.b2base_thresh       = b2base_thresh; 
E_hit_cal.num_E2_valid        = num_E2_valid; 
E_hit_cal.num_E1_valid        = num_E1_valid; 
E_hit_cal.num_hits_b2base     = num_hits_b2base; 
E_hit_cal.num_hits_no_b2base  = num_hits_no_b2base; 
E_hit_cal.reward_prob_per_frame   = reward_prob_per_frame; 


%hit_data
E_hit_data                    = E_hit_cal;
E_hit_data.T_hits             = T_hits; 
E_hit_data.E2_valid           = E2_valid; 
E_hit_data.E1_valid           = E1_valid; 
E_hit_data.T_idxs_no_b2base   = T_idxs_no_b2base;
E_hit_data.T_idxs_b2base      = T_idxs_b2base;    

end

function [score] = ...
        score_calibration(reward_prob, reward_per_frame_range)
    
score = 0; 
within_range = (reward_prob >= reward_per_frame_range(1)) && (reward_prob <= reward_per_frame_range(2)); 
if(~within_range)
    score = min(abs(reward_prob - reward_per_frame_range)); 
end

end

function [score, E2_score, E1_score] = score_calibration_2target(E1_reward_prob, E2_reward_prob, reward_per_frame_range)
[E2_score] = ...
    score_calibration(E2_reward_prob, reward_per_frame_range);
[E1_score] = ...
    score_calibration(E1_reward_prob, reward_per_frame_range);
score = E1_score + E2_score; 
end

function plot_cal_hits(hit_data,...
    E1_mean_analyze, E2_mean_analyze, ...
    cursor_obs, title_str, save_path)

cursor_amp = (max(cursor_obs)-min(cursor_obs));
cursor_offset = cursor_amp/10; 
max_cursor = max(cursor_obs); 


c1 = hit_data.T_hits; 
c2 = hit_data.E2_valid; 
c3 = hit_data.E1_valid;
valid_hit_idxs = hit_data.T_idxs_b2base;
num_valid_hits = hit_data.num_hits_b2base;
% hit_data.T_hits             = T_hits; 
% hit_data.E2_valid           = E2_valid; 
% hit_data.E1_valid           = E1_valid; 
% hit_data.T_idxs_no_b2base   = T_idxs_no_b2base;
% hit_data.T_idxs_b2base      = T_idxs_b2base;


h =figure; hold on;
scatter(c1, ones(length(c1),1)*max_cursor + cursor_offset, 15, 'r'); %plot(cursor_obs-cursor_offset, 'k'); 
scatter(c2, ones(length(c2),1)*max_cursor + 2*cursor_offset, 15, 'g'); %plot(cursor_obs-cursor_offset, 'k'); 
scatter(c3, ones(length(c3),1)*max_cursor + 3*cursor_offset, 15, 'b'); %plot(cursor_obs-cursor_offset, 'k'); 
plot(cursor_obs); 
hline(hit_data.T); 
plot(E1_mean_analyze-cursor_amp); 
plot(E2_mean_analyze-2*cursor_amp); 
xlabel('frame'); 
title([title_str ' hits with b2base: ' num2str(num_valid_hits)]); 
legend({'c1', 'c2 - E1 cond', 'c3 - E2 cond', 'cursor', 'E1 mean', 'E2 subord mean'}); 
if(num_valid_hits > 0)
    vline(valid_hit_idxs); 
end
saveas(h, save_path);

end
