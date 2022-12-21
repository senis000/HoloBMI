function simulate_all_bmi (Experiments, rows)
% Experiments is a table imported from the parquet file 

variable_n_data_name = 'file_training';
variable_target_info_name = 'file_target_info';
variable_baseline = 'file_baseline';

folder_re_sim = 'C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\re_simulations_BMI';
folder_sim = 'C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\simulations_target2_BMI';
var_n_data = {};
var_target_info = {};
var_baseline = {};
NVI12 = nan(size(rows));
NVI13 = nan(size(rows));
NVI16 = nan(size(rows));
NVI17 = nan(size(rows));
NVI20 = nan(size(rows));
NVI22 = nan(size(rows));
did_it_worked = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);
number_hits_t1 = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);
number_hits_t2 = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);
cursor_occupancy = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);
cursor_occupancy_t2 = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);
cursor_occupancy_baseline = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);
cursor_occupancy_t2_baseline = table(rows, NVI12, NVI13, NVI16, NVI17, NVI20, NVI22);

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
        tic
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
                mkdir(folder_to_save);
            end
            if ~ exist(folder_to_save_t2_cal)
                mkdir(folder_to_save_t2_cal);
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
            len_data = data.frame;
            clear data;
            %% obtain back2base of baseline
            n_f_file = fullfile(folder_raw, file_baseline);
            [num_hits_no_b2base_t1, len_base] = obtain_back2base_baseline(n_f_file, target_info.E1_base, ...
                target_info.E2_base, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, ...
                f0_init_slide, target_info.T1, target_info.E1_thresh, target_info.E2_subord_thresh);
            cursor_occupancy_baseline(row, animal{1}) = {num_hits_no_b2base_t1/len_base*len_data};
            close all
            %% simulate BMI
            [~, num_simulated_hits, ~, occupancy_cursor] = ...
                sim_bmi_vE1strict_fb_corrected(n_data, ...
                task_settings, target_info, folder_to_save);
            close all
            if number_hits == num_simulated_hits
                did_it_worked(row, animal{1}) = {1};
            else
                did_it_worked(row, animal{1}) = {0};
            end
            number_hits_t1(row, animal{1}) = {number_hits};
            cursor_occupancy(row, animal{1}) = {occupancy_cursor};
            %% simulation otherway around
            % obtain with baseline the calibration 
            A_file = fullfile(folder_raw, 'roi_data.mat');
            onacid_bool = 0;
            E1_base = target_info.E2_base;
            E2_base = target_info.E1_base;
            [target_info_t2_path, ~, ~, num_hits_no_b2base_t2] = ...
                baseline2target_fb_objective(n_f_file, A_file, onacid_bool,  ...
                E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
                prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, folder_to_save_t2_cal, ...
                cursor_zscore_bool, f0_init_slide, E2mE1_prctile, fb_settings);
            target_info_t2 = load(target_info_t2_path); 
            close all
            [~, num_simulated_t2_hits, ~, occupancy_cursor_t2] = ...
                sim_bmi_vE1strict_fb_corrected(n_data, ...
                task_settings, target_info_t2, folder_to_save_t2);
            number_hits_t2(row, animal{1}) = {num_simulated_t2_hits};
            cursor_occupancy_t2(row, animal{1}) = {occupancy_cursor_t2};
            cursor_occupancy_t2_baseline(row, animal{1}) = {num_hits_no_b2base_t2/len_base*len_data};
            close all
        end
        disp(['Finishing animal: ', animal, ' day: ', int2str(rows(row)), ' lasted: ', num2str(toc)]);
    end
    
    
end

parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\sim_worked.parquet', did_it_worked)
parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\cursor_occupancy.parquet', cursor_occupancy)
parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\cursor_occupancy_t2.parquet', cursor_occupancy_t2)
parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\hits_t2.parquet', number_hits_t2)
parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\hits.parquet', number_hits_t1)
parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\cursor_occupancy_baseline.parquet', cursor_occupancy_baseline)
parquetwrite('C:\Users\Nuria\Documents\DATA\holoBMI\curated_data\cursor_occupancy_t2_baseline.parquet', cursor_occupancy_t2_baseline)

end
