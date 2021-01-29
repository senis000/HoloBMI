function [num_hits_no_b2base] = obtain_back2base_baseline(n_f_file, ...
    E1_base, E2_base, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, f0_init_slide, T, E1_thresh, ...
    E2_subord_thresh)


%1) Select Temporal: BMI E data from baseline 

load(n_f_file); 
f_base = baseActivity; 
%n_f_file - contains matrix, neural fluorescence from baseline file, num_samples X num_neurons_baseline 
f_base(:,isnan(f_base(1,:))) = [];
f_base = f_base.'; %num_samples X num_neurons

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
E2_sel = E_id==2; 
E2_sel_idxs = find(E2_sel); 


%%
%4) Decoder information
%Decoder information

[decoder, E1_proj, E2_proj, E1_norm, E2_norm] = def_decoder(num_neurons, E1_sel, E2_sel);

%%
%First process f0: 
if(f0_win_bool)    
    %Calculate f0 as in BMI: 
    
    if f0_init_slide
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
    else
        num_samples = size(f_raw,1);
        f0 = zeros(num_samples-f0_win+1, num_neurons); 
        f0(1,:) = mean(f_raw(1:f0_win, :), 1);
        for i = 2:length(f0)
            f0(i,:) = f0(i-1,:)*((f0_win-1)/f0_win) + f_raw((i+f0_win-1), :)/f0_win; 
        end
        %Truncate data based on the f0_win:
        f_postf0 = f_raw(f0_win:end, :); 
    end
else
    f_postf0 = f_raw; 
    f0 = f0_mean; 
end

%%
%Second, smooth f:
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


%%
%Third, compute dff and dff_z:
dff = (f_smooth-f0)./f0;


%%
n_analyze = dff;
 
valid_idxs  = find(~isnan(n_analyze(:,1)));
n_analyze   = n_analyze(valid_idxs, :); 
analyze_cov = cov(n_analyze);
analyze_mean = nanmean(n_analyze); %takes mean along dim 1.  n_analyze is num_samples X num_neurons.

%%
%Inputs: 
%frames_per_reward_range
%cov_bool
E2_subord_mean          = zeros(num_E2,1);
E2_subord_std           = zeros(num_E2,1); 
E1_analyze              = n_analyze(:,E1_sel); 
E2_analyze              = n_analyze(:,E2_sel); 
for E2_i = 1:num_E2
    subord_sel                      = E2_sel;
    subord_sel(E2_sel_idxs(E2_i))   = 0; 
    E2_subord_mean(E2_i)            = mean(analyze_mean(subord_sel));     
    var_i                           = subord_sel'*analyze_cov*subord_sel; 
    E2_subord_std(E2_i)             = sqrt(var_i);     
end

E2_sum_analyze = sum(E2_analyze,2); 

%signals needed for target detection:
cursor_obs                      = n_analyze*decoder; 
E1_mean_analyze                 = mean(E1_analyze,2); 
[E2_dom_samples, E2_dom_sel]    = max(E2_analyze, [], 2);
E2_subord_mean_analyze          = (E2_sum_analyze - E2_dom_samples)/(num_E2-1);

%
% h = figure;
% hist(mean(E1_analyze,2)); 
% vline(E1_mean); 
% vline(E1_mean+E1_std); 

%%
%1) E2-E1 > alpha
c1 = find(cursor_obs >= T); 
%2) E1 < mu
c2 = find(E1_mean_analyze <= E1_thresh);
%3) E2_subord > mu (anded with previous constraint)
%For each idx, subtract the 
c3 = find(E2_subord_mean_analyze >= E2_subord_thresh(E2_dom_sel)); 
hit_idxs_no_b2base = intersect(intersect(c1, c2), c3);

num_hits_no_b2base = length(hit_idxs_no_b2base);
end


function [decoder, E1_proj, E2_proj, E1_norm, E2_norm] = ...
    def_decoder(num_neurons, E1_sel, E2_sel)

E1_proj = zeros(num_neurons, 1); 
E1_proj(E1_sel) = 1;
E1_norm = sum(E1_sel); %can replace with vector norm.  
% disp('E1 proj'); 
E1_proj = E1_proj/E1_norm;
% E1_proj

E2_proj = zeros(num_neurons, 1); 
E2_proj(E2_sel) = 1; 
E2_norm = sum(E2_sel); 
% disp('E2 proj')
E2_proj = E2_proj/E2_norm;
% E2_proj

% disp('decoder:')
decoder = E2_proj - E1_proj;
end

 