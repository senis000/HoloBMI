%results_110419.m
%%
%export_fig
%hline
addpath('/Users/vivekathalye/Dropbox/Code/MatlabExchange/hline_vline'); 
addpath('/Users/vivekathalye/Dropbox/Code/export_fig/export_fig_11.19.15'); 

%%
save_dir = '/Users/vivekathalye/Dropbox/Data/holo_nov2019/fb_performance_check_110419/plots';

%%
data_home   = '/Users/vivekathalye/Dropbox/Data/holo_nov2019/fb_performance_check_110419'
dates       = {'191030', '191031', '191101', '191102', '191103', '191104'}; 
num_days    = length(dates); 

stim_idxs   = [2 4 6]; 
num_stim    = length(stim_idxs); 
nostim_idxs = [1 3 5]; 
num_nostim  = length(nostim_dates); 

animals = {'NVI12', 'NVI13', 'NVI16'}; 
num_animals = length(animals); 

base_mat = zeros(num_animals, num_days); 
bmi_mat = zeros(num_animals, num_days); 
%Want matrices of size: num_animals x num_dates: 
%base_stim
%base_nostim
%bmi_stim
%bmi_nostim

% base_cell   = cell(num_animals, num_days); 
bmi_cell    = cell(num_animals, num_days); 


%%
%STIM DATA: 
for i_day = 1:num_days
    disp('day: '); 
    disp(i_day); 
    for i_animal = 1:num_animals
        disp('animal: '); 
        disp(i_animal); 
        
        dir_i = fullfile(data_home, dates{i_day}, animals{i_animal}); 
%         exist(dir_i)

        dir_files = dir(dir_i);

        num_base = 0; 
        num_bmi = 0; 
        base_paths = {};
        bmi_paths = {};
        for i = 1:length(dir_files)
            if ~isempty(strfind(dir_files(i).name, 'target_calibration_ALL'))
                num_base = num_base + 1; 
                base_paths{num_base}    = fullfile(dir_i, dir_files(i).name); 
            elseif ~isempty(strfind(dir_files(i).name, 'BMI_online'))
                num_bmi = num_bmi + 1; 
                bmi_paths{num_bmi}      = fullfile(dir_i, dir_files(i).name); 
            end
        end

        base_path   = base_paths{end}; 
        bmi_path    = bmi_paths{end}; 

        base_data   = load(base_path); 
        bmi_data    = load(bmi_path); 
        
        %ASSIGN:
        base_mat(i_animal, i_day) = base_data.num_valid_hits; 
        bmi_mat(i_animal, i_day) = bmi_data.data.selfTargetCounter;
        
        bmi_cell{i_animal, i_day} = ...
            bmi_data.data.selfHits; 
    end
end

%%
base_mat_scaled = base_mat*40/15; 

%%
animal_colors = {'r', 'g', 'k'}; 
h = figure;
hold on;
for i = 1:num_animals
%     plot(base_mat_scaled(i,:), animal_colors{i}); 
    plot(bmi_mat(i,:), '.-', 'MarkerSize', 15, 'color', animal_colors{i}); 
end

set(gca,'TickDir','out');
export_fig(h, fullfile(save_dir, 'perf_across_days.eps')); 


%%
base_stim   = base_mat_scaled(:,stim_idxs); 
bmi_stim    = bmi_mat(:,stim_idxs); 
base_nostim = base_mat_scaled(:,nostim_idxs); 
bmi_nostim  = bmi_mat(:,nostim_idxs); 

%%
bmi_stim_mean = mean(bmi_stim, 2)
bmi_nostim_mean = mean(bmi_nostim, 2)


bar_c = [1 2 3]*3; 
bar_width = 3
xdelta_bar = 0.5; 

h = figure;
hold on; 
for i_animal =1:num_animals
    %no stim: 
    bar(bar_c(i_animal)-xdelta_bar, bmi_nostim_mean(i_animal), 'k'); 
    x_data = ones(num_days/2,1)*bar_c(i_animal)-xdelta_bar
    scatter(x_data, bmi_nostim(i_animal, :), 70, 'b', 'filled'); 
    
    bar(bar_c(i_animal)+xdelta_bar, bmi_stim_mean(i_animal), 'r'); 
    
    x_data = ones(num_days/2,1)*bar_c(i_animal)+xdelta_bar
    scatter(x_data, bmi_stim(i_animal, :), 70, 'b', 'filled');
end

hline(40/15*7); 

ylabel('number of hits'); 
set(gca, 'XTickLabel', []); 
set(gca,'TickDir','out');
xlabel('animal')

title('black: pretrain with random reward, red: pretrain with paired stim and reward'); 
export_fig(h, fullfile(save_dir, 'stim_vs_nostim_animal.eps')); 
%%
bmi_stim_vec = bmi_stim(:)
bmi_nostim_vec = bmi_nostim(:)

bmi_stim_mean = mean(bmi_stim_vec); 
bmi_nostim_mean = mean(bmi_nostim_vec); 

bar_width = 0.7; 
xdelta_bar = 0.5; 
h = figure;
hold on;

bar(0-xdelta_bar, bmi_nostim_mean, bar_width, 'k')
x_data = ones(length(bmi_stim_vec),1)*0 -xdelta_bar
scatter(x_data, bmi_nostim_vec, 70, 'b', 'filled'); 
    
bar(0+xdelta_bar, bmi_stim_mean, bar_width, 'r')
x_data = ones(length(bmi_stim_vec),1)*0 +xdelta_bar
scatter(x_data, bmi_stim_vec, 70, 'b', 'filled'); 
ylabel('num hits'); 
hline(40/15*7); 

[p,h] = ranksum(bmi_nostim_vec, bmi_stim_vec)

title(['black: pretrain random reward, red: pretrain paired stim reward, ranksum p: ' num2str(p)]); 

set(gca,'TickDir','out');

export_fig(h, fullfile(save_dir, 'stim_vs_nostim_pool.eps')); 

%%
%Loop over and do: 


%%
%Do an example: 
selfHits = bmi_cell{1,1}; 

num_valid = length(~isnan(selfHits))

%%
h = figure; 
plot(selfHits); 

%%
num_samples         = 10000; 
smooth_filt         = ones(num_samples,1)/num_samples; 
selfHits_smooth     = conv(selfHits, smooth_filt, 'valid'); 

%%
h = figure;
plot(selfHits_smooth)

    
%%
num_samples         = 10000; 
smooth_filt         = ones(num_samples,1)/num_samples; 

smooth_hits_cell = cell(num_animals, num_days); 
for i_day = 1:num_days
    for i_animal = 1:num_animals
        selfHits = bmi_cell{i_animal,i_day}; 
        selfHits_smooth = conv(selfHits, smooth_filt, 'valid');
        
        smooth_hits_cell{i_animal, i_day} = selfHits_smooth; 
    end
end

%%
stim_r      = zeros(size(selfHits_smooth));
nostim_r    = zeros(size(selfHits_smooth));

for i_day = 1:num_days
    for i_animal = 1:num_animals
        if sum(i_day == stim_idxs) > 0
            stim_r = stim_r + ...
                smooth_hits_cell{i_animal, i_day}; 
        else
            nostim_r = nostim_r + ...
                smooth_hits_cell{i_animal, i_day};             
        end
        
    end
end
stim_r = stim_r/(num_animals*3); 
nostim_r = nostim_r/(num_animals*3); 

stim_r_trunc    = stim_r(1:70000)*30*60; 
nostim_r_trunc  = nostim_r(1:70000)*30*60; 

x_vec = (1:70000)/(30*60); 

%%
h = figure;
plot(smooth_hits_cell{1,2})
%%

h = figure;
hold on;
plot(x_vec, stim_r_trunc, 'r', 'LineWidth', 2); 
plot(x_vec, nostim_r_trunc, 'k', 'LineWidth', 2); 
legend({'stim reward pretrain', 'random reward pretrain'}); 
xlabel('time (min)'); 
ylabel('hits per min'); 
set(gca,'TickDir','out');
title('within session learning, each curve is average of 9 animal-sessions'); 

export_fig(h, fullfile(save_dir, 'within_session_learning.eps'));
export_fig(h, fullfile(save_dir, 'within_session_learning.png'));



%%
h = figure; 
hold on;
for i_day = 1:num_days
    for i_animal = 1:num_animals
        if sum(i_day == stim_idxs) > 0
            plot(x_vec, smooth_hits_cell{i_animal, i_day}(1:70000)*30*60, 'r');
        end
%         else
%             plot(smooth_hits_cell{i_animal, i_day}, 'k'); 
%         end        
    end
end
plot(x_vec, stim_r_trunc, 'r', 'LineWidth', 4); 

%%
h = figure; 
hold on;
for i_day = 1:num_days
    for i_animal = 1:num_animals
        if sum(i_day == stim_idxs) > 0
%             plot(x_vec, smooth_hits_cell{i_animal, i_day}(1:70000)*30*60, 'r');
        
        else
            plot(smooth_hits_cell{i_animal, i_day}(1:70000)*30*60, 'k'); 
        end        
    end
end
plot(x_vec, stim_r_trunc, 'k', 'LineWidth', 4); 
