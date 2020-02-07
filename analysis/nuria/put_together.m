function put_together(basemin, toplot, ELtime)

%% nargins and other inputs
if nargin <1
    basemin = 5;
end

if nargin <2
    toplot = true;
end
if nargin <3
    ELtime = 14;
end

%preparing the structure
dirpath = 'C:/Data/Nuria/HoloBMI/experiments/';
dirpathanalysis = 'C:/Data/Nuria/HoloBMI/analysis/';

%% create the vars
animal(1).name = 'NVI12';
animal(2).name = 'NVI13';
animal(3).name = 'NVI16';
animal(4).name = 'NVI17';
animal(5).name = 'NVI20';
animal(6).name = 'NVI22';
%just to bring the values of the days when we did the experiments
animal(1).days.holobmi = [190930, 191001, 191005, 191007, 191009, 191011, 191013, 191025]; 
animal(2).days.holobmi = [190930, 191001, 191003, 191005, 191007, 191009, 191011, 191013, 191025]; 
animal(3).days.holobmi = [190930, 191001, 191005, 191007, 191009, 191011, 191025];

animal(1).days.E3 = [191004, 191026, 191028];
animal(2).days.E3 = [191004, 191026, 191028];
animal(3).days.E3 = [191004, 191026, 191028];

animal(1).days.rr = [191006, 191010, 191014];
animal(2).days.rr = [191006, 191010, 191014];
animal(3).days.rr = [191006, 191010, 191014, 191024];

animal(1).days.E2holo = [191008, 191012, 191027];
animal(2).days.E2holo = [191008, 191012, 191027];
animal(3).days.E2holo = [191008, 191012, 191027];

animal(1).days.frr = [191030, 191101, 191103];
animal(2).days.frr = [191030, 191101, 191103];
animal(3).days.frr = [191030, 191101, 191103];

animal(1).days.fhbmi = [191031, 191102, 191104];
animal(2).days.fhbmi = [191031, 191102, 191104];
animal(3).days.fhbmi = [191031, 191102, 191104];


animal(4).days.holobmi = [191106, 191108, 191110, 191112, 191114, 191116, 191120, 191122];
animal(5).days.holobmi = [191106, 191108, 191110, 191112, 191114, 191116, 191118, 191120];
animal(6).days.holobmi = [191106, 191108, 191110, 191114, 191116, 191118, 191120, 191122];

animal(4).days.E3 = [191109, 191115, 191121];
animal(5).days.E3 = [191109, 191115, 191121];
animal(6).days.E3 = [191109, 191115, 191121];

animal(4).days.rr = [191111, 191117, 191123];
animal(5).days.rr = [191105, 191117, 191122];
animal(6).days.rr = [191105, 191111, 191117];

animal(4).days.E2holo = [191107, 191113, 191119];
animal(5).days.E2holo = [191107, 191113, 191119];
animal(6).days.E2holo = [191107, 191113, 191119];

animal(4).days.frr = [191124, 191126, 191128];
animal(5).days.frr = [191124, 191126, 191128];
animal(6).days.frr = [191124, 191126, 191128];

animal(4).days.fhbmi = [191125, 191127, 191129];
animal(5).days.fhbmi = [191125, 191127, 191130];
animal(6).days.fhbmi = [191125, 191127, 191129];

%% MEASURES OF LEARNING definitions by animal
%measures of learning
ML.aa.TH = nan(length(animal),6); %total hits 
ML.aa.THbase = nan(length(animal),6); %total hits in baseline
ML.aa.THg = nan(length(animal),6); %(total hits -baseline hits) / baseline hits
ML.aa.HPM = nan(length(animal),6); %hits per minute
ML.aa.HPMg = nan(length(animal),6); %(hits per minute -baseline hitspm) / baseline hitspm

%measures of learning when baseline is first 3 min of BMI
ML.aa.THbg = nan(length(animal),6); %total hits - baselinehits / baseline hits
ML.aa.HPMbg = nan(length(animal),6); %hits per minute - hitsbaseline / hitsbaseline

%measures of learning based on gain from early/late
ML.aa.HitGain = nan(length(animal),6); % ((hits Early -hits Late )/ hits early)
ML.aa.HPMGain = nan(length(animal),6); % ((hpm early - hpm late )/ hpm early)

%measures based on cursor occupancy
ML.aa.CO = nan(length(animal),6); % % of occupancy
ML.aa.COPer = nan(length(animal),6); % % of occupancy
ML.aa.COGain = nan(length(animal),6); % gain in occupancy from E1-E2
ML.aa.COHits = nan(length(animal),6); % Gain of hits from E2 -> E1

%hits per min curve
ML.aa.HPMcurve = nan(length(animal),6, 39); %hits per minute

%% Measures of learning by session
% calculate how many sessions
totalDays = zeros(length(animal),1);
for aa=1:length(animal)
    typeDays = fieldnames(animal(aa).days);
    for tt=1:length(typeDays)
        totalDays(aa) = totalDays(aa) + length(getfield(animal(aa).days, typeDays{tt}));
    end
        
end
totalSessions = max(totalDays);

%measures of learning
ML.ss.TH = nan(totalSessions,6); %total hits 
ML.ss.THbase = nan(length(animal),6); %total hits in baseline
ML.ss.THg = nan(totalSessions,6); %(total hits -baseline hits) / baseline hits
ML.ss.HPM = nan(totalSessions,6); %hits per minute
ML.ss.HPMg = nan(totalSessions,6); %(hits per minute -baseline hitspm) / baseline hitspm

%measures of learning when baseline is first 3 min of BMI
ML.ss.THbg = nan(totalSessions,6); %total hits - baselinehits / baseline hits
ML.ss.HPMbg = nan(totalSessions,6); %hits per minute - hitsbaseline / hitsbaseline

%measures of learning based on gain from early/late
ML.ss.HitGain = nan(totalSessions,6); % ((hits Early -hits Late )/ hits early)
ML.ss.HPMGain = nan(totalSessions,6); % ((hpm early - hpm late )/ hpm early)

%measures based on cursor occupancy
ML.ss.CO = nan(totalSessions,6); % % of occupancy
ML.ss.COPer = nan(totalSessions,6); % % of occupancy
ML.ss.COGain = nan(totalSessions,6); % gain in occupancy from E1-E2
ML.ss.COHits = nan(totalSessions,6); % Gain of hits from E2 -> E1

%hits per min curve
ML.ss.HPMcurve = nan(totalSessions,6, 39); %hits per minute

%% 

%% Obtain measures
session = zeros(1,6); %number of sessions per each experiment
for aa=1:length(animal)
    typeDays = fieldnames(animal(aa).days);
    
    for tt=1:length(typeDays)
        daystoGet = animal(aa).days.(typeDays{tt});
        aux.TH = nan(1,length(daystoGet));
        aux.THbase = nan(1,length(daystoGet));
        aux.THg = nan(1,length(daystoGet));
        aux.HPM = nan(1,length(daystoGet));
        aux.HPMg = nan(1,length(daystoGet));
        aux.THbg = nan(1,length(daystoGet));
        aux.HPMbg = nan(1,length(daystoGet));
        aux.HitGain = nan(1,length(daystoGet));
        aux.HPMGain = nan(1,length(daystoGet));
        aux.CO = nan(1,length(daystoGet));
        aux.COPer = nan(1,length(daystoGet));
        aux.COGain = nan(1,length(daystoGet));
        aux.COHits = nan(1,length(daystoGet));
        aux.hpm_curve = nan(length(daystoGet), 39);
        for dd = 1:length(daystoGet)
            session(tt) = session(tt) + 1;
            %this is all to find the bmi_online files
            folder.dirPath = fullfile(dirpath,num2str(daystoGet(dd)),animal(aa).name);
            folder.simPath = fullfile(folder.dirPath,'Simulations','BMIsim');
            folder.simPathcont = fullfile(folder.dirPath,'Simulations','BMIcont');
            if ~exist(folder.simPath, 'dir')
                mkdir(folder.simPath)
            end
            if ~exist(folder.simPathcont, 'dir')
                mkdir(folder.simPathcont)
            end
            filesPath = dir(folder.dirPath);
            ibmi=1;
            clear('bmiFiles');
            for ff=3:length(filesPath)
                if strcmp(filesPath(ff).name(1:5),'BMI_o')
                    bmiFiles{ibmi} = filesPath(ff).name;
                    ibmi = ibmi+1;
                end
                if strcmp(filesPath(ff).name(1:3),'wor')
                    workfile = filesPath(ff).name;
                end
            end
            %if there is less than 2 bmi_online something is wrong
            if length(bmiFiles)<2
                disp('STHAAAAP Something is wrong with this experiment')
                disp(folder.dirPath)
                return
            else
                fileBMIOpen = fullfile(folder.dirPath, bmiFiles{end});
                fileWorksOpen = fullfile(folder.dirPath, workfile);
                load(fileWorksOpen, 'task_settings','target_info_path','target_cal_ALL_path');
                calibrat = load(fileWorksOpen, 'n_f_file', 'A_file', 'frames_per_reward_range');
                load(fileBMIOpen, 'data');
                [~,aux.bname,aux.bext] = fileparts(target_cal_ALL_path);
                [~,aux.tiname,aux.tiext] = fileparts(target_info_path);
                new_tc_path = fullfile(folder.dirPath, strcat(aux.bname, aux.bext));
                new_ti_path = fullfile(folder.dirPath, strcat(aux.tiname, aux.tiext));
                baseData = load(new_tc_path, 'cursor_obs', 'hits', 'hits_valid');
                target_info = load(new_ti_path);
            end
            % experimental variables
            binsExp = linspace(3600, data.frame, 40);  % for bins starting after the 2 min
            %because if I start at 0 I'm biasing the results to have lower
            %hits/min during the first 5min of baseline I start from task_settings.prefix_win
            binsExpb = linspace(task_settings.prefix_win, data.frame+task_settings.prefix_win, 42); % for bins 

starting at 0 ~ task_settings.prefix_win
      
            basehits = nansum(baseData.hits_valid)*40/15; % hits if baseline would be 40min
            basehpm = nansum(baseData.hits_valid)/15; %hpm of baseline
            
            % calculate/store measures
            aux.TH(dd) = nansum(data.selfHits);
            aux.THbase(dd) = nansum(baseData.hits_valid);
            aux.THg(dd) = (nansum(data.selfHits) - basehits)/basehits*100;
            aux.hpm_curve(dd,:) = histcounts(find(data.selfHits==1), binsExp);
            aux.HPM(dd) = nanmean(aux.hpm_curve(dd,:));
            aux.HPMg(dd) = (aux.HPM(dd) - basehpm) / basehpm*100;
            
            % to calculate bmi from the very beginning
            task_settings.calibration.f0_init_slide = 1;
            % create or retrieve the simulation of cursor during BMI
            % if the plots dir doesn't exist most likely the sim file
            % doesn't exist either
            if ~exist(fullfile(folder.simPath,'plots'),'dir')
                [sim_saved_path] = sim_bmi_vE1strict_fb(data.bmiAct, task_settings, target_info, 

folder.simPath);
                close('all')
            else
                isim = 1;
                filesSimPath = dir(folder.simPath);
                for ff=3:length(filesSimPath)
                    if strcmp(filesSimPath(ff).name(1:3),'sim')
                        simFiles{isim} = filesSimPath(ff).name;
                        isim = isim+1;
                    end
                end
                sim_saved_path = fullfile(folder.simPath, simFiles{end});
            end
            simData = load(sim_saved_path, 'cursor_obs', 'valid_hit_idxs');
            
            % measures with baseline define at the begining of BMI
            aux_hpmb_curve = histcounts(simData.valid_hit_idxs, binsExpb);
            if nansum(aux_hpmb_curve(1:basemin)) > 0
                basehitsb = nansum(aux_hpmb_curve(1:basemin))*(41-basemin)/basemin; %total hits if baseline 

would be 42 - basemin minutes
                basehpmb = nanmean(aux_hpmb_curve(1:basemin)); %hpm 
                aux.THbg(dd) = (nansum(aux_hpmb_curve(basemin+1:end)) - basehitsb)/ basehitsb*100;
                aux.HPMbg(dd) = (nanmean(aux_hpmb_curve(basemin+1:end)) - basehpmb)/ basehpmb*100;

                % measures of early vs late
                aux.HitGain(dd) = (nansum(aux_hpmb_curve(end-ELtime+1:end)) - nansum(aux_hpmb_curve

(1:ELtime)))/ nansum(aux_hpmb_curve(1:ELtime))*100;
                aux.HPMGain(dd) = (nanmean(aux_hpmb_curve(end-ELtime+1:end)) - nanmean(aux_hpmb_curve

(1:ELtime)))/ nanmean(aux_hpmb_curve(1:ELtime))*100;
            end
            %lets deal with the cursor occupancy now
            task_settings.calibration.f0_init_slide = 0;
            if ~exist(fullfile(folder.simPathcont,'plots'),'dir')
                %create the contrary calibration
                calibrat.E2mE1_prctile = 98; 
                if ~exist('fb_settings','var')
                    [fb_settings] = define_fb_audio_settings();
                end
                % obtain the names of the files 
                [~,auxbasename,auxext] = fileparts(calibrat.n_f_file);
                [~,auxroiname,auxroiext] = fileparts(calibrat.A_file);
                calibrat.basename = fullfile(folder.dirPath, strcat(auxbasename, auxext));
                calibrat.roiname = fullfile(folder.dirPath, strcat(auxroiname, auxroiext));
                % calibrate with the oposite E1 and E2
                [target_rev_path, target_rev_ALL_path, ~] = baseline2target_vE1strict_fb(calibrat.basename, 

calibrat.roiname, 0,  ...
                    target_info.E2_base, target_info.E1_base, calibrat.frames_per_reward_range, 0, ...
                    task_settings.prefix_win, 1, task_settings.f0_win, 1, task_settings.dff_win, 

folder.simPathcont, ...
                    0, task_settings.calibration.f0_init_slide, calibrat.E2mE1_prctile, fb_settings);
                close('all')
                target_rev_info = load(target_rev_path);
                task_rev_settings = load(target_rev_ALL_path);
                % simulate with the oposite E1 and E2
                task_rev_settings.calibration.f0_win_bool = task_rev_settings.f0_win_bool;
                task_rev_settings.calibration.f0_init_slide = task_rev_settings.f0_init_slide;
                task_rev_settings.back2BaseFrameThresh = task_rev_settings.back2BaseFramesThresh;
                save(fullfile(folder.simPathcont,'task_rev.mat'), 'task_rev_settings');
                [sim_saved_path_cursor] = sim_bmi_vE1strict_fb(data.bmiAct, task_rev_settings, target_rev_info, 

folder.simPathcont);
                close('all')
            else
                isimco = 1;
                filesSimPathcont = dir(folder.simPathcont);
                for ff=3:length(filesSimPathcont)
                    if strcmp(filesSimPathcont(ff).name(1:3),'sim')
                        simFilescont{isimco} = filesSimPathcont(ff).name;
                        isimco = isimco+1;
                    end
                end
                sim_saved_path_cursor = fullfile(folder.simPathcont, simFilescont{end});
                load(fullfile(folder.simPathcont,'task_rev.mat'), 'task_rev_settings');
            end
            simDataCursor = load(sim_saved_path_cursor, 'cursor_obs', 'valid_hit_idxs');
            length_rev_cursor = length(find(simDataCursor.cursor_obs > task_rev_settings.T));
            length_cursor = length(find(data.cursor > target_info.T1));
            % measure with the simulated cursor of opposite E2
            aux.CO(dd) = length_cursor/data.frame;
            aux.COPer(dd) = length_cursor/(length_rev_cursor+length_cursor);
            aux.COHits(dd) = nansum(data.selfHits)/(nansum(data.selfHits)+length

(simDataCursor.valid_hit_idxs));
            aux.COGain(dd) = aux.COHits(dd)/(nansum(baseData.hits_valid)/(nansum(baseData.hits_valid) + 

task_rev_settings.num_valid_hits ));
            if toplot
                h = figure();
                plot(1.5:41.5, aux_hpmb_curve)
                hold on
                plot(3.5:41.5, aux.hpm_curve(dd,:))
                legend({'simulated', 'original'}); 
                xlabel('time (min)'); 
                ylabel('Hits'); 
                title('Hits per min');
                im_path = fullfile(folder.simPath, 'plots','hpm_test.png'); 
                saveas(h, im_path);
                close('all')

                h2 = figure();
                minc = min(min(baseData.cursor_obs), min(data.cursor));
                maxc = max(max(baseData.cursor_obs), max(data.cursor));
                binsed = minc:0.05:maxc;
                hisb = histcounts(baseData.cursor_obs, binsed)/sum(~isnan(baseData.cursor_obs));
                hise = histcounts(data.cursor, binsed)/sum(~isnan(data.cursor));
                ba = bar(binsed(2:end),hisb);
                hold on;
                bb =bar(binsed(2:end),hise);
                alpha(ba,0.5);
                alpha(bb,0.5);
                vline(target_info.T1, 'r')
                vline(-task_rev_settings.T, 'k')
                legend('Baseline','BMI');
                dirplot = fullfile(folder.dirPath, 'plots');
                if ~exist(dirplot,'dir')
                    mkdir(dirplot)
                end
                im_path2 = fullfile(dirplot,'base_exp_cursor.png'); 
                saveas(h2, im_path2);
                im_path3 = fullfile(dirpathanalysis, 'plots', strcat(typeDays{tt}, animal(aa).name, '_', 

num2str(dd), '.png')); 
                saveas(h2, im_path3);
                close('all')
                    
            end
            %update values for session
            ML.ss.TH(session(tt),tt) = aux.TH(dd); %total hits 
            ML.ss.THbase(session(tt),tt) = aux.THbase(dd); %total hits in baseline
            ML.ss.THg(session(tt),tt) = aux.THg(dd); %(total hits -baseline hits) / baseline hits
            ML.ss.HPM(session(tt),tt) = aux.HPM(dd); %hits per minute
            ML.ss.HPMg(session(tt),tt) = aux.HPMg(dd); %(hits per minute -baseline hitspm) / baseline hitspm

            %measures of learning when baseline is first 3 min of BMI
            ML.ss.THbg(session(tt),tt) = aux.THbg(dd); %total hits - baselinehits / baseline hits
            ML.ss.HPMbg(session(tt),tt) = aux.HPMbg(dd); %hits per minute - hitsbaseline / hitsbaseline

            %measures of learning based on gain from early/late
            ML.ss.HitGain(session(tt),tt) = aux.HitGain(dd); % ((hits Early -hits Late )/ hits early)
            ML.ss.HPMGain(session(tt),tt) = aux.HPMGain(dd); % ((hpm early - hpm late )/ hpm early)

            %measures based on cursor occupancy
            ML.ss.CO(session(tt),tt) =  aux.CO(dd);  % % of occupancy
            ML.ss.COPer(session(tt),tt) =  aux.COPer(dd);  % % of occupancy
            ML.ss.COGain(session(tt),tt) = aux.COGain(dd); % gain in occupancy from E1-E2
            ML.ss.COHits(session(tt),tt) = aux.COHits(dd);  % Gain of hits from E2 -> E1

            %hits per min curve
            ML.ss.HPMcurve(session(tt),tt,:) = aux.hpm_curve(dd,:); %hits per minute
        end
        ML.aa.TH(aa,tt) = nanmean(aux.TH); %total hits 
        ML.aa.THbase(aa,tt) = nanmean(aux.THbase); %total hits in baseline
        ML.aa.THg(aa,tt) = nanmean(aux.THg); %(total hits -baseline hits) / baseline hits
        ML.aa.HPM(aa,tt) = nanmean(aux.HPM); %hits per minute
        ML.aa.HPMg(aa,tt) = nanmean(aux.HPMg); %(hits per minute -baseline hitspm) / baseline hitspm

        %measures of learning when baseline is first 3 min of BMI
        ML.aa.THbg(aa,tt) = nanmean(aux.THbg); %total hits - baselinehits / baseline hits
        ML.aa.HPMbg(aa,tt) = nanmean(aux.HPMbg); %hits per minute - hitsbaseline / hitsbaseline

        %measures of learning based on gain from early/late
        ML.aa.HitGain(aa,tt) = nanmean(aux.HitGain); % ((hits Early -hits Late )/ hits early)
        ML.aa.HPMGain(aa,tt) = nanmean(aux.HPMGain); % ((hpm early - hpm late )/ hpm early)

        %measures based on cursor occupancy
        ML.aa.CO(aa,tt) = nanmean(aux.CO); % % of occupancy
        ML.aa.COPer(aa,tt) = nanmean(aux.COPer); % % of occupancy
        ML.aa.COGain(aa,tt) = nanmean(aux.COGain); % gain in occupancy from E1-E2
        ML.aa.COHits(aa,tt) = nanmean(aux.COHits); % Gain of hits from E2 -> E1

        %hits per min curve
        ML.aa.HPMcurve(aa,tt,:) = nanmean(aux.hpm_curve,1); %hits per minute

    end
end
MLpathanalysis = fullfile(dirpathanalysis, 'ML.mat');
save(MLpathanalysis, 'ML', 'animal')

%% plot results and save
bar(nanmean(ML.aa.TH))

%% plots
subplot(221)
bar(nanmean(ML.aa.TH(:,1:4),1),'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.TH(:,1:4),1)/sqrt(length(6));
errorbar(nanmean(ML.aa.TH(:,1:4),1), SEM, '.k')
title('Total Hits')
ylabel('Hits')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest2(ML.aa.TH(:,1),ML.aa.TH(:,2));
[~,p3] = ttest2(ML.aa.TH(:,1),ML.aa.TH(:,3));
[~,p4] = ttest2(ML.aa.TH(:,1),ML.aa.TH(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([0,50])


subplot(222)
bar(nanmean(ML.aa.HPM(:,1:4),1),'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.HPM(:,1:4),1)/sqrt(length(6));
errorbar(nanmean(ML.aa.HPM(:,1:4),1), SEM, '.k')
title('Hits per min')
ylabel('Hpm')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest2(ML.aa.HPM(:,1),ML.aa.HPM(:,2));
[~,p3] = ttest2(ML.aa.HPM(:,1),ML.aa.HPM(:,3));
[~,p4] = ttest2(ML.aa.HPM(:,1),ML.aa.HPM(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([0,1.3])

subplot(223)
bar(nanmean(ML.aa.THg(:,1:4),1),'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.THg(:,1:4),1)/sqrt(length(6));
errorbar(nanmean(ML.aa.THg(:,1:4),1), SEM, '.k')
title('Hits gain')
ylabel('Increase Hits (%)')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest2(ML.aa.THg(:,1),ML.aa.THg(:,2));
[~,p3] = ttest2(ML.aa.THg(:,1),ML.aa.THg(:,3));
[~,p4] = ttest2(ML.aa.THg(:,1),ML.aa.THg(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([-70,140])


subplot(224)
bar(nanmean(ML.aa.HPMg(:,1:4),1),'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.HPMg(:,1:4),1)/sqrt(length(6));
errorbar(nanmean(ML.aa.HPMg(:,1:4),1), SEM, '.k')
title('HPM Gain')
ylabel('Increase Hpm (%)')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest(ML.aa.HPMg(:,1),ML.aa.HPMg(:,2));
[~,p3] = ttest(ML.aa.HPMg(:,1),ML.aa.HPMg(:,3));
[~,p4] = ttest(ML.aa.HPMg(:,1),ML.aa.HPMg(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([-70,140])

%% individual mice
figure()
for i=1:6
    subplot(3,2,i)
    bar(ML.aa.HPM(i,1:4),'FaceColor', uint8([110 110 110]));
    xlabel(animal(i).name)
    ylabel('HPM')
    ylim([0,1])
    xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
end

%% learning from baseline

figure
bar([nanmean(ML.aa.THbase(:,1)/15,1), nanmean(ML.aa.TH(:,1)/40,1)],'FaceColor', uint8([110 110 110]));
SEM = [nanstd(ML.aa.THbase(:,1)/15,1), nanstd(ML.aa.TH(:,1)/40,1)]/sqrt(length(6));
hold on
errorbar([nanmean(ML.aa.THbase(:,1)/15,1), nanmean(ML.aa.TH(:,1)/40,1)], SEM, '.k')
title('Hits per min')
ylabel('Hits per min')
xticklabels({'Baseline', 'BMI'})
[~,p2] = ttest(ML.aa.THbase(:,1)/15,ML.aa.TH(:,1)/40);
sigstar([1,2],p2);
ylim([0,1.3])

%% plots of CO
figure()
subplot(221)
bar(nanmean(ML.aa.CO(:,1:4),1)*100,'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.CO(:,1:4),1)*100/sqrt(length(6));
errorbar(nanmean(ML.aa.CO(:,1:4),1)*100, SEM, '.k')
title('Ocupancy T1')
ylabel('Occupancy (%)')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest(ML.aa.CO(:,1),ML.aa.CO(:,2));
[~,p3] = ttest(ML.aa.CO(:,1),ML.aa.CO(:,3));
[~,p4] = ttest(ML.aa.CO(:,1),ML.aa.CO(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([0,4])


subplot(222)
bar(nanmean(ML.aa.COGain(:,1:4),1),'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.COGain(:,1:4),1)/sqrt(length(6));
errorbar(nanmean(ML.aa.COGain(:,1:4),1), SEM, '.k')
title('T1/T2 hit ratio')
ylabel('ratio')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest2(ML.aa.COGain(:,1),ML.aa.COGain(:,2));
[~,p3] = ttest2(ML.aa.COGain(:,1),ML.aa.COGain(:,3));
[~,p4] = ttest2(ML.aa.COGain(:,1),ML.aa.COGain(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([0,2.7])

subplot(223)
bar(nanmean(ML.aa.COHits(:,1:4),1),'FaceColor', uint8([110 110 110]));
hold on;
SEM = nanstd(ML.aa.COHits(:,1:4),1)/sqrt(length(6));
errorbar(nanmean(ML.aa.COHits(:,1:4),1), SEM, '.k')
title('T1/T2 increase')
ylabel('ratio')
xticklabels({'Stim', 'E3', 'NR', 'RR', 'FRR', 'FS'})
[~,p2] = ttest2(ML.aa.COHits(:,1),ML.aa.COHits(:,2));
[~,p3] = ttest2(ML.aa.COHits(:,1),ML.aa.COHits(:,3));
[~,p4] = ttest2(ML.aa.COHits(:,1),ML.aa.COHits(:,4));

sigstar([1,2],p2);
sigstar([1,3],p3);
sigstar([1,4],p4);
ylim([0,1.3])
