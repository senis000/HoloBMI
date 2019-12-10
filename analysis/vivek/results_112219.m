%results_112219

%%
%export_fig
%hline
addpath('/Users/vivekathalye/Dropbox/Code/MatlabExchange/hline_vline'); 
addpath('/Users/vivekathalye/Dropbox/Code/export_fig/export_fig_11.19.15');

%%
get_data_from_server = 0;
server_dir = ...
    '/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1'
%     '/Volumes/costa-labshare/Vivek/VivekNuria/HoloBMI_2ndround';
save_dir = ...
    '/Users/vivekathalye/Dropbox/Data/bmi_test_nov2019/fb_performance_check_112219'
plot_dir = ...
    fullfile(save_dir, 'plots'); 
mkdir(plot_dir); 

dates       = {'2019-11-19', '2019-11-21', '2019-11-22'}; 
%days        = {'D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', ...
    %'D09', 'D10', 'D11', 'D12', 'D13', 'D14'};
num_days    = length(dates); 

% animals     = {'NVI17', 'NVI20', 'NVI21', 'NVI22'}; 
animals     = {'NY127'}; 
num_animals     = length(animals); 


%Need two files: 
%baseline calibration file
%BMI online file
%if can't find either one, then that date is invalid.  
%
if get_data_from_server
	%loop date, animal, day
    valid_over_ad   = ones(num_animals, num_days); 
    for i_day = 1:num_days
        date_i  = dates{i_day};
%         day_i   = days{i_day}; 
        for i_animal = 1:num_animals
            animal_i    = animals{i_animal}; 
            date_day_animal_dir = fullfile(animal_i, date_i); 
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
                if ~isempty(strfind(dir_files(i).name, 'BMI_cal_ALL'))
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
%Manual checks: 

%%
cal_path = ...
    fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NVI18', '2019-11-22', 'BMI_cal_ALL_20191122T100403.mat'); 
cal_test = load(cal_path); 
data_path = fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NVI18', '2019-11-22', 'BMI_online191122T104206.mat'); 
test = load(data_path)

%%
cal_path = ...
    fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-22', 'BMI_cal_ALL_20191122T112510.mat'); 
cal_test = load(cal_path); 
data_path = fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-22', 'BMI_online191122T120838.mat'); 
test = load(data_path)

%%
% cal_path = ...
%     fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-21', 'BMI_online191121T102019.mat'); 
% cal_test = load(cal_path); 
data_path = fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-21', 'BMI_online191121T102019.mat'); 
test = load(data_path)

%%
% cal_path = ...
%     fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-21', 'BMI_online191121T102019.mat'); 
% cal_test = load(cal_path); 
data_path = fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-21', 'BMI_online191121T104607.mat'); 
test = load(data_path)

%%
data_path = fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-21', 'BMI_online191121T111214.mat'); 
test = load(data_path)

%%
data_path = fullfile('/Volumes/costa-locker/Ines/2photon/Training/ISO+BMI 1', 'NY127', '2019-11-19', 'BMI_online191119T102404.mat'); 
test = load(data_path)

%%
E2_valid        = find(~isnan(test.data.selfHits)); 
num_E2_valid    = length(E2_valid); 
t_E2 = (1:num_E2_valid)/1800; 

E1_valid        = find(~isnan(test.data.E1Hits)); 
num_E1_valid    = length(E1_valid)
t_E1 = (1:num_E1_valid)/1800; 

E2hit_cumsum = cumsum(test.data.selfHits(E2_valid)); 
E1hit_cumsum = cumsum(test.data.E1Hits(E1_valid)); 
h = figure; 
hold on;
plot(t_E1, E1hit_cumsum); 
plot(t_E2, E2hit_cumsum); 
legend('unrewarded E1 hits', 'rewarded E2 hits'); 
% legend('unrewarded E2 hits (high freq)', 'rewarded E1 hits (low freq)'); 
% legend('unrewarded E1 hits (low freq)', 'rewarded E2 hits (high freq)'); 

% cal_num_E1 = cal_test.cal.target.E1_hit_cal.num_hits_b2base; 
% cal_num_E2 = cal_test.cal.target.E2_hit_cal.num_hits_b2base; 
% hline(cal_num_E1*40/15); 

xlabel('minutes'); 
ylabel('number of hits'); 
