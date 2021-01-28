function [did_it_worked] = simulate_all_bmi (Experiments, rows)
% Experiments is a table imported from the parquet file 

variable_n_data_name = 'file_training';
variable_target_info_name = 'file_target_info';
variable_baseline = 'file_baseline';

folder_sim = 'C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\re_simulations_BMI';
var_n_data = {};
var_target_info = {};
var_baseline = {};
did_it_worked = zeros(size(rows), 'logical');
% variables for simulation
% baseline_frameRate = num_base_samples/(15*60);
% sec_per_reward_range = [120 85]; 
% frames_per_reward_range = sec_per_reward_range*baseline_frameRate;
% E2mE1_prctile = 98; 
% target_on_cov_bool = 0;
% prefix_win = 40;
% f0_win_bool = 1;
% f0_win = 2*60*ceil(frameRate);
% dff_win_bool = 1;
% dff_win = 4;
% cursor_zscore_bool = 0;
% f0_init_slide = 0;
% [fb_settings]   = define_fb_audio_settings();


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
            folder_to_save = fullfile(folder_sim, int2str(rows(row)), animal, day);
            if ~ exist(folder_to_save)
                mkdir(folder_to_save)
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
            clear data;
            %% simulate BMI
            [~, num_simulated_hits, valid_hit_idxs] = ...
                sim_bmi_vE1strict_fb_corrected(n_data, ...
                task_settings, target_info, folder_to_save);
            if number_hits == num_simulated_hits
                did_it_worked(row) = true;
            end
            close all
            %% simulation otherway around
            % obtain with baseline the calibration 
    %         n_f_file = fullfile(folder_raw, file_baseline);
    %         A_file = fullfile(folder_raw, 'roi_data.mat');
    %         onacid_bool = 0;
    %         E1_base = target_info.E2_base;
    %         E2_base = target_info.E1_base;
    %         [target_info_path, target_cal_ALL_path, fb_cal] = ...
    %             baseline2target_fb_objective(n_f_file, A_file, onacid_bool,  ...
    %             E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
    %             prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
    %             cursor_zscore_bool, f0_init_slide, E2mE1_prctile, fb_settings);
        
        end
    end
    
    
end
