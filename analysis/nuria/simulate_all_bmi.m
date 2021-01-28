function [did_it_worked, number_hits_t2] = simulate_all_bmi (Experiments, rows)
% Experiments is a table imported from the parquet file 

variable_n_data_name = 'file_training';
variable_target_info_name = 'file_target_info';
variable_baseline = 'file_baseline';

folder_re_sim = 'C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\re_simulations_BMI';
folder_sim = 'C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\simulations_target2_BMI';
var_n_data = {};
var_target_info = {};
var_baseline = {};
did_it_worked = zeros(size(rows), 'logical');
difference_hits = zeros(size(rows), 'double');
number_hits_t2 = zeros(size(rows), 'single');
% variables for simulation

sec_per_reward_range = [120 85]; 
baseline_frameRate = 30;
frameRate = 30;
frames_per_reward_range = sec_per_reward_range*baseline_frameRate;
E2mE1_prctile = 98; 
target_on_cov_bool = 0;
prefix_win = 40;
f0_win_bool = 1;
f0_win = 2*60*ceil(frameRate);
dff_win_bool = 1;
dff_win = 4;
cursor_zscore_bool = 0;
f0_init_slide = 0;
[fb_settings]   = define_fb_audio_settings();


for str_var = Experiments.Properties.VariableNames
    if size(strfind(str_var{1}, variable_n_data_name),1)
        var_n_data{end+1} = str_var{1};
    end
    if size(strfind(str_var{1}, variable_target_info_name),1)
        var_target_info{end+1} = str_var{1};
    end
    if size(strfind(str_var{1}, variable_baseline),1)
        var_baseline{end+1} = str_var{1};
    end
end

number_rows = size(Experiments, 1);
for row=1:number_rows
    for animal={('NVI12'), ('NVI13'), ('NVI16'), ('NVI17'), ('NVI20'), ('NVI22')}
        %% open files
        aux_day = [Experiments(row, 'x__Day___1_round__').Variables, ...
            Experiments(row, 'x__Day___2_round__').Variables];
        day = aux_day(~ismissing(aux_day));
        file_data = NaN;
        for str_var = var_n_data
            if size(strfind(str_var{1}, animal{1}),1)
                file_data = Experiments(row, str_var{1}).Variables;
            end
        end
        if ~ ismissing(file_data)
            for str_var = var_target_info
                if size(strfind(str_var{1}, animal{1}),1)
                    file_target_info = Experiments(row, str_var{1}).Variables;
                end
            end
            for str_var = var_baseline
                if size(strfind(str_var{1}, animal{1}),1)
                    file_baseline = Experiments(row, str_var{1}).Variables;
                end
            end
        
            folder_raw = fullfile ("C:\Users\Nuria\Documents\DATA\holoBMI\raw_data\", ...
                int2str(rows(row)), animal, day);
            folder_to_save = fullfile(folder_re_sim, int2str(rows(row)), animal, day);
            folder_to_save_t2 = fullfile(folder_sim, int2str(rows(row)), animal, day);
            folder_to_save_t2_cal = fullfile(folder_sim, int2str(rows(row)), animal, day, 'calibration');
            if ~ exist(folder_to_save)
                mkdir(folder_to_save)
            end
            if ~ exist(folder_to_save_t2_cal)
                mkdir(folder_to_save_t2_cal)
            end
            load(fullfile(folder_raw, file_data), 'data');
            try
                load(fullfile(folder_raw, 'workspace.mat'), 'task_settings');
            catch
                disp('using previous workspace')
            end
            target_info = load(fullfile(folder_raw, file_target_info));       
            n_data = data.bmiAct;
            number_hits = data.selfTargetVTACounter;
            location_hits = find(data.selfVTA)';
            clear data;
            %% simulate BMI
            [~, num_simulated_hits, valid_hit_idxs] = ...
                sim_bmi_vE1strict_fb_corrected(n_data, ...
                task_settings, target_info, folder_to_save);
            if number_hits == num_simulated_hits
                did_it_worked(row) = true;
                difference_hits(row) = sum(location_hits - valid_hit_idxs);
            end
            close all
            %% simulation otherway around
            % obtain with baseline the calibration 
            n_f_file = fullfile(folder_raw, file_baseline);
            A_file = fullfile(folder_raw, 'roi_data.mat');
            onacid_bool = 0;
            E1_base = target_info.E2_base;
            E2_base = target_info.E1_base;
            [target_info_t2_path, ~, ~] = ...
                baseline2target_fb_objective(n_f_file, A_file, onacid_bool,  ...
                E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
                prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, folder_to_save_t2_cal, ...
                cursor_zscore_bool, f0_init_slide, E2mE1_prctile, fb_settings);
            target_info_t2 = load(target_info_t2_path); 
            close all
            [~, num_simulated_t2_hits, ~] = ...
                sim_bmi_vE1strict_fb_corrected(n_data, ...
                task_settings, target_info_t2, folder_to_save_t2);
            number_hits_t2(row) = num_simulated_t2_hits;
            close all
        end
    end
    
    
end
