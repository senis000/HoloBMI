%results_111319.m
%%
%export_fig
%hline
addpath('/Users/vivekathalye/Dropbox/Code/MatlabExchange/hline_vline'); 
addpath('/Users/vivekathalye/Dropbox/Code/export_fig/export_fig_11.19.15'); 

%%
get_data_from_server = 0;
server_dir = ...
    '/Volumes/costa-labshare/Vivek/VivekNuria/HoloBMI_2ndround';
save_dir = ...
    '/Users/vivekathalye/Dropbox/Data/holo_nov2019/fb_performance_check_111319'
plot_dir = ...
    fullfile(save_dir, 'plots'); 
mkdir(plot_dir); 

dates       = {'191105', '191106', '191107', '191108', '191109', '191110', '191111', '191112'}; 
days        = {'D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08'};
num_days    = length(days); 

animals     = {'NVI17', 'NVI20', 'NVI21', 'NVI22'}; 
num_animals     = length(animals); 


%Need two files: 
%baseline calibration file
%BMI online file
%if can't find either one, then that date is invalid.  
%%
if get_data_from_server
	%loop date, animal, day
    valid_over_ad   = ones(num_animals, num_days); 
    for i_day = 1:num_days
        date_i  = dates{i_day}; 
        day_i   = days{i_day}; 
        for i_animal = 1:num_animals
            animal_i    = animals{i_animal}; 
            date_day_animal_dir = fullfile(date_i, animal_i, day_i); 
            dir_i    = fullfile(server_dir, date_day_animal_dir); 
            
            %Get the files: 
            dir_files = dir(dir_i)
            num_base = 0; 
            num_bmi = 0; 
            base_paths  = {};
            base_files  = {};
            bmi_paths   = {};
            bmi_files   = {};
            
            for i = 1:length(dir_files)
                if ~isempty(strfind(dir_files(i).name, 'target_calibration_ALL'))
                    num_base = num_base + 1;
                    
                    base_files{num_base}    = dir_files(i).name; 
                    base_paths{num_base}    = fullfile(dir_i, dir_files(i).name); 
                elseif ~isempty(strfind(dir_files(i).name, 'BMI_online'))
                    num_bmi = num_bmi + 1;
                    
                    bmi_files{num_bmi}      = dir_files(i).name; 
                    bmi_paths{num_bmi}      = fullfile(dir_i, dir_files(i).name); 
                end
            end
            %Make save directory: 
            save_dir_i = fullfile(save_dir, date_day_animal_dir); 
            mkdir(save_dir_i);            
            if(num_base ==0 || num_bmi ==0)
                valid_over_ad(i_animal, i_day) = 0; 
            else
                if(num_base > 0)
                    base_file   = base_files{end}; 
                    base_path   = base_paths{end};
                    %Save the data in save directory:  
                    copyfile(base_path, save_dir_i);
                end
                if(num_bmi > 0)
                    bmi_file    = bmi_files{end}; 
                    bmi_path    = bmi_paths{end}; 
                    %Save the data in save directory:                      
                    copyfile(bmi_path, save_dir_i);                    
                end             
            end
        end
    end
end
%%
base_mat        = zeros(num_animals, num_days); 
bmi_mat         = zeros(num_animals, num_days); 
%Want matrices of size: num_animals x num_dates: 
% base_cell   = cell(num_animals, num_days); 
bmi_cell        = cell(num_animals, num_days); 
valid_over_ad   = zeros(num_animals, num_days); 

for i_day = 1:num_days
    date_i  = dates{i_day}; 
    day_i   = days{i_day}; 
    for i_animal = 1:num_animals
        animal_i    = animals{i_animal}; 
        date_day_animal_dir = fullfile(date_i, animal_i, day_i); 
        dir_i    = fullfile(save_dir, date_day_animal_dir); 

        %Get the files: 
        dir_files = dir(dir_i)
        num_base = 0; 
        num_bmi = 0; 
        base_paths  = {};
        base_files  = {};
        bmi_paths   = {};
        bmi_files   = {};

        for i = 1:length(dir_files)
            if ~isempty(strfind(dir_files(i).name, 'target_calibration_ALL'))
                num_base = num_base + 1;

                base_files{num_base}    = dir_files(i).name; 
                base_paths{num_base}    = fullfile(dir_i, dir_files(i).name); 
            elseif ~isempty(strfind(dir_files(i).name, 'BMI_online'))
                num_bmi = num_bmi + 1;

                bmi_files{num_bmi}      = dir_files(i).name; 
                bmi_paths{num_bmi}      = fullfile(dir_i, dir_files(i).name); 
            end
        end
          
        if(num_base > 0  && num_bmi > 0)
            valid_over_ad(i_animal, i_day) = 1;             
        end
        
        if(num_base > 0)
            base_file   = base_files{end}; 
            base_path   = base_paths{end};
            base_data   = load(base_path); 

            base_mat(i_animal, i_day) = base_data.num_valid_hits; 
        end
        if(num_bmi > 0)
            bmi_file    = bmi_files{end}; 
            bmi_path    = bmi_paths{end}; 
            bmi_data    = load(bmi_path); 

            bmi_mat(i_animal, i_day) = bmi_data.data.selfTargetCounter;
            bmi_cell{i_animal, i_day} = ...
                bmi_data.data.selfHits; 
        end
    end
end

%%
bmi_mat_backup = bmi_mat
%%
%Filter out invalid data: 
valid_over_ad(2,7) = 0; 
% bmi_mat = bmi_mat.*valid_over_ad;
bmi_mat(find(~valid_over_ad)) = NaN

%%
base_mat_scaled = base_mat*40/15; 


%%
animal_colors = {'r', 'g', 'b', 'k'}; 
h = figure;
hold on;
for i = 1:num_animals
%     plot(base_mat_scaled(i,:), animal_colors{i}); 
    plot(bmi_mat(i,:), '.-', 'MarkerSize', 15, 'color', animal_colors{i}); 
end
legend(animals)
%%
stim_idxs       = [2 4 6 8]
nostim_idxs     = [1 7]

stim_valid      = valid_over_ad(:,stim_idxs); 
nostim_valid    = valid_over_ad(:,nostim_idxs); 

base_stim   = base_mat_scaled(:,stim_idxs); 
bmi_stim    = bmi_mat(:,stim_idxs); 
base_nostim = base_mat_scaled(:,nostim_idxs); 
bmi_nostim  = bmi_mat(:,nostim_idxs); 

%%
bmi_nostim

%%
bmi_nostim(find(nostim_valid))
%%
bmi_stim_mean       = nanmean(bmi_stim, 2)
bmi_nostim_mean     = nanmean(bmi_nostim, 2)


bar_c = [1 2 3 4]*3; 
bar_width = 3
xdelta_bar = 0.5; 

h = figure;
hold on; 
for i_animal =1:num_animals

    bar(bar_c(i_animal)-xdelta_bar, bmi_nostim_mean(i_animal), 'k'); 
    y_data = bmi_nostim(i_animal, :)
    x_data = ones(length(y_data),1)*bar_c(i_animal)-xdelta_bar
    scatter(x_data, y_data, 70, 'b', 'filled'); 
    
    
    bar(bar_c(i_animal)+xdelta_bar, bmi_stim_mean(i_animal), 'r'); 
    y_data = bmi_stim(i_animal, :)
    x_data = ones(length(y_data),1)*bar_c(i_animal)+xdelta_bar
    scatter(x_data, bmi_stim(i_animal, :), 70, 'b', 'filled');
end

hline(40/15*7); 

ylabel('number of hits'); 
set(gca, 'XTickLabel', []); 
set(gca,'TickDir','out');
xlabel('animal')

title('NVI 17 20 21 22 black: pretrain with random reward, red: pretrain with paired stim and reward'); 
export_fig(h, fullfile(plot_dir, 'stim_vs_nostim_animal.eps')); 
export_fig(h, fullfile(plot_dir, 'stim_vs_nostim_animal.png')); 

%%
%Remove animal 3
bmi_stim_remove     = bmi_stim([1 2 4], :)
bmi_nostim_remove   = bmi_nostim([1 2 4], :)

bmi_stim_vec = bmi_stim_remove(:)
bmi_nostim_vec = bmi_nostim_remove(:)

bmi_stim_mean = nanmean(bmi_stim_vec); 
bmi_nostim_mean = nanmean(bmi_nostim_vec); 

bar_width = 0.7; 
xdelta_bar = 0.5; 
h = figure;
hold on;

bar(0-xdelta_bar, bmi_nostim_mean, bar_width, 'k')
x_data = ones(length(bmi_nostim_vec),1)*0 -xdelta_bar
scatter(x_data, bmi_nostim_vec, 70, 'b', 'filled'); 
    
bar(0+xdelta_bar, bmi_stim_mean, bar_width, 'r')
x_data = ones(length(bmi_stim_vec),1)*0 +xdelta_bar
scatter(x_data, bmi_stim_vec, 70, 'b', 'filled'); 
ylabel('num hits'); 
hline(40/15*7); 

[p,h] = ranksum(bmi_nostim_vec, bmi_stim_vec)

title(['black: pretrain random reward, red: pretrain paired stim reward, ranksum p: ' num2str(p)]); 

set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, 'stim_vs_nostim_pool_remove.eps')); 
export_fig(h, fullfile(plot_dir, 'stim_vs_nostim_pool_remove.png')); 

%%

bmi_stim_vec = bmi_stim(:)
bmi_nostim_vec = bmi_nostim(:)

bmi_stim_mean = nanmean(bmi_stim_vec); 
bmi_nostim_mean = nanmean(bmi_nostim_vec); 

bar_width = 0.7; 
xdelta_bar = 0.5; 
h = figure;
hold on;

bar(0-xdelta_bar, bmi_nostim_mean, bar_width, 'k')
x_data = ones(length(bmi_nostim_vec),1)*0 -xdelta_bar
scatter(x_data, bmi_nostim_vec, 70, 'b', 'filled'); 
    
bar(0+xdelta_bar, bmi_stim_mean, bar_width, 'r')
x_data = ones(length(bmi_stim_vec),1)*0 +xdelta_bar
scatter(x_data, bmi_stim_vec, 70, 'b', 'filled'); 
ylabel('num hits'); 
hline(40/15*7); 

[p,h] = ranksum(bmi_nostim_vec, bmi_stim_vec)

title(['black: pretrain random reward, red: pretrain paired stim reward, ranksum p: ' num2str(p)]); 

set(gca,'TickDir','out');

export_fig(h, fullfile(plot_dir, 'stim_vs_nostim_pool.eps')); 
export_fig(h, fullfile(plot_dir, 'stim_vs_nostim_pool.png')); 
