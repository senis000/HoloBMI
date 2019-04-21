function [h, early_late_mean, early_late_sem, early_vs_late_pval] = ...
    plot_early_vs_late(data_str, data_over_var_cond, valid_over_var_cond, colors_over_cond, ...
    early_idxs, late_idxs, ranksum_side, ...
    ylim_bool, ylim_vec, save_dir)
%INPUT:
%ranksum_side is a string, one of these: {'left', 'right', ''}.  
%'' means ranksum will not be a sided test
%'left' means test late > early
%'right' means test early > late
%
%STATS:
% mu early, mu late, ranksum test
% num consistent targets

num_var         = size(data_over_var_cond, 1);
num_cond        = size(data_over_var_cond, 2);

%select early
sel_idxs = early_idxs;
[early_pool,early_data, early_sel] = ...
    sel_data(data_over_var_cond, valid_over_var_cond, sel_idxs, num_var, num_cond);
disp('early')
early_data
early_pool

%select late
sel_idxs = late_idxs;
[late_pool, late_data, late_sel] = ...
    sel_data(data_over_var_cond, valid_over_var_cond, sel_idxs, num_var, num_cond);

disp('late')
late_data
late_pool

% if(mean(early_pool) < mean(late_pool))
%     ranksum_side    = 'left';
% else
%     ranksum_side    = 'right';
% end

%1) Plot each condition in Gray
%2) Plot each cond in Colors

% %GRAY
% colors = 0.7*ones(num_cond, 3);
% [h, ...
%     early_vs_late_pval, early_late_mean, early_late_sem, ...
%     num_consistent_cond] = util_plot_early_vs_late(early_pool, early_data, early_idxs, ...
%     late_pool, late_data, late_idxs, colors);
% export_fig(h, fullfile(save_dir, [data_str ' early_late_over_cond_gray' '.eps']));

% %COLORS
colors = colors_over_cond;
[h, ...
    early_vs_late_pval, early_late_mean, early_late_sem, ...
    num_consistent_cond] = util_plot_early_vs_late(...
    valid_over_var_cond, ...
    early_pool, early_data, early_idxs, ...
    late_pool, late_data, late_idxs, ...
    colors, ranksum_side, ylim_bool, ylim_vec);
export_fig(h, fullfile(save_dir, [data_str ' early_late_over_cond' '.eps']));

end

function [data_sel_pool,data_sel_mat, sel] = ...
    sel_data(data_over_var_cond, valid_over_var_cond, sel_idxs, num_var, num_cond)
    sel_mat             = zeros(num_var, num_cond);
    sel_mat(sel_idxs, 1:num_cond) = 1;
    sel                 = valid_over_var_cond & sel_mat
    data_sel_pool       = data_over_var_cond(sel);
    data_sel_mat        = data_over_var_cond(sel_idxs,:);
end

% for animal_idx = 1:num_animals
%     plot([1.1 1.9], [early_mean late_mean], '.-', 'LineWidth', 1.5, 'color', 0.6*[1 1 1], 'Marker', 'o'); 
% end
% errbar([1 2], [early_mean_pool late_mean_pool], [early_sem_pool late_sem_pool],  'LineWidth', 2, 'Color', 'k');
% plot([1 2], [early_mean_pool late_mean_pool], 'LineWidth', 2, 'Color', 'k'); 

function [h, ...
    early_vs_late_pval, early_late_mean, early_late_sem, ...
    num_consistent_cond] = util_plot_early_vs_late(...
    valid_over_var_cond, ...
    early_pool, early_data, early_idxs, ...
    late_pool, late_data, late_idxs, colors_over_cond, ...
    ranksum_side, ylim_bool, ylim_vec);
    
    num_cond = size(valid_over_var_cond, 2);     

    %Early vs Late Pval 
    [early_vs_late_pval,H] = ranksum(early_pool, late_pool, 'tail', ranksum_side);       
    %ranksum_side - comes from earlier analysis

    early_pool
    late_pool
    [early_sem, early_mean]  = sem(early_pool(:)', ones(size(early_pool(:)')))
    [late_sem, late_mean]    = sem(late_pool(:)', ones(size(late_pool(:)')))
    early_late_mean = [early_mean late_mean];
    early_late_sem  = [early_sem late_sem];

    % function [sem_result, mean_result] = sem(data_over_var_cond, valid_over_var_cond)
    % %DESCRIPTION:
    % %Calculates standard error of mean: s/sqrt(n)
    % %INPUT:
    % %data - num_variables X num observations
    % %valid - num_variables X num_observations


    h = figure;
    hold on;
    %--------------------------------------------------------------------------
    %PLOT EACH COND
    x_data = [1.1 1.9];
    num_consistent_cond = 0;
    for cond_idx = 1:num_cond
        early_valid_i       = find(valid_over_var_cond(early_idxs, cond_idx));
        early_data_i        = early_data(early_valid_i, cond_idx);

        late_valid_i    = find(valid_over_var_cond(late_idxs, cond_idx));
        late_data_i     = late_data(late_valid_i, cond_idx);

        mu_early_i  = mean(early_data_i);
        mu_late_i   = mean(late_data_i);

        y_data = [mu_early_i mu_late_i];
        plot(x_data, y_data, '.-', 'LineWidth', 1.5, 'color', colors_over_cond(cond_idx,:), 'Marker', 'o');

        if(strcmp(ranksum_side, 'left'))
            if(mu_early_i <= mu_late_i)
                num_consistent_cond = ...
                    num_consistent_cond + 1;
            end
        else
            if(mu_early_i >= mu_late_i)
                num_consistent_cond = ...
                    num_consistent_cond + 1;
            end        
        end
    end
    
    %--------------------------------------------------------------------------
    %PLOT MEAN+SEM    
    x_data = [1 2];
    errbar(1:2,early_late_mean, early_late_sem, 'LineWidth', 2, 'color', [0 0 0]);
    plot(x_data, early_late_mean, '.-', 'color', 'k', 'LineWidth', 2);
    xlim([0 3]);
    
    
    if(ylim_bool)
        ylim(ylim_vec);
    end
    
    set(gca, 'XTick', 1:2);
    set(gca, 'XTickLabel', {'early', 'late'});
    set(gca, 'TickDir', 'out');
    set(gca, 'box', 'off');
    varargout = sigstar({[1 2]}, [early_vs_late_pval]);

    title_str = ['early: ' num2str(early_mean) ', late: ' num2str(late_mean) ...
        ' ranksum: ' num2str(early_vs_late_pval) ' num cond consistent: ' num2str(num_consistent_cond)];
    title(title_str);    
end
    

% function [h, early_late_mean, early_late_sem, early_vs_late_pval] = ...
%     plot_early_vs_late(data_str, data_over_target_epoch, valid_over_target_epoch, save_dir)
% %STATS:
% % mu early, mu late, ranksum test
% % num consistent targets
% 
% %%
% % Debug: 
% % data_over_target_epoch = fa_result.pop_var_over_target_epoch.shared_main;
% % data_over_target_epoch = fa_result.pop_var_over_target_epoch.shared;
% %%
% num_epochs      = size(data_over_target_epoch, 2);
% NUM_TARGETS     = size(data_over_target_epoch, 1);
% 
% num_epochs_compare = 2;
% early_epochs    = 1:1+num_epochs_compare;
% late_epochs     = (num_epochs-num_epochs_compare):num_epochs;
% 
% %Select valid early epochs
% early_sel   = valid_over_target_epoch;
% early_sel(:, early_epochs(end)+1:end) = 0;
% early_sel   = early_sel == 1; %convert to logical
% %Select valid late epochs
% late_sel    = valid_over_target_epoch;
% late_sel(:, 1:(late_epochs(1)-1)) = 0;
% late_sel   = late_sel == 1; %convert to logical
% 
% %Select data
% early_data  = data_over_target_epoch(:,1:1+num_epochs_compare);
% late_data   = data_over_target_epoch(:,(end-num_epochs_compare):end);
% 
% early_pool  = data_over_target_epoch(early_sel);
% late_pool   = data_over_target_epoch(late_sel); 
% 
% % early_data  = data_over_target_epoch(:,1:1+num_epochs_compare);
% % late_data   = data_over_target_epoch(:,(end-num_epochs_compare):end);
% % early_pool  = early_data(:);
% % late_pool   = late_data(:);
% 
% %--------------------------------------------------------------------------
% 
% if(mean(early_pool) < mean(late_pool))
%     ranksum_side    = 'left';
% else
%     ranksum_side    = 'right';
% end
% 
% %Early vs Late Pval 
% [early_vs_late_pval,H] = ranksum(early_pool, late_pool, 'tail', ranksum_side);       
% %ranksum_side - comes from earlier analysis
% 
% [early_sem early_mean]  = sem(early_pool', ones(size(early_pool')))
% [late_sem late_mean]    = sem(late_pool', ones(size(late_pool')))
% early_late_mean = [early_mean late_mean];
% early_late_sem  = [early_sem late_sem];
% 
% % function [sem_result, mean_result] = sem(data, valid);
% % %function sem_result = sem(data);
% % %DESCRIPTION:
% % %Calculates standard error of mean: s/sqrt(n)
% % %INPUT:
% % %data - num_variables X num observations
% % %valid - num_variables X num_observations
% 
% %--------------------------------------------------------------------------
% %PLOT MEAN+SEM
% x_data = [1 2];
% h = figure;
% hold on;
% errbar(1:2,early_late_mean, early_late_sem, 'LineWidth', 2, 'color', [0 0 0]);
% plot(x_data, early_late_mean, '.-', 'color', 'k', 'LineWidth', 2);
% xlim([0 3]);
% 
% %--------------------------------------------------------------------------
% %PLOT EACH TARGET
% colors = vivek_colors();
% x_data = [1.1 1.9];
% num_consistent_targets = 0;
% for target_idx = 1:NUM_TARGETS
%     early_valid_i     = find(valid_over_target_epoch(target_idx, early_epochs));
%     early_data_i    = early_data(target_idx, early_valid_i);
%     %early_data_i    = early_data(target_idx, :);
%     
%     late_valid_i    = find(valid_over_target_epoch(target_idx, late_epochs));
%     late_data_i     = late_data(target_idx, late_valid_i);
%     %late_data_i     = late_data(target_idx, :);
%     
%     mu_early_i  = mean(early_data_i);
%     mu_late_i   = mean(late_data_i);
%     
%     y_data = [mu_early_i mu_late_i];
%     plot(x_data, y_data, 'color', colors(target_idx,:), 'LineWidth', 1.25);
%     
%     if(strcmp(ranksum_side, 'left'))
%         if(mu_early_i <= mu_late_i)
%             num_consistent_targets = ...
%                 num_consistent_targets + 1;
%         end
%     else
%         if(mu_early_i >= mu_late_i)
%             num_consistent_targets = ...
%                 num_consistent_targets + 1;
%         end        
%     end
% end
% 
% set(gca, 'XTick', 1:2);
% set(gca, 'XTickLabel', {'early', 'late'});
% set(gca, 'TickDir', 'out');
% set(gca, 'box', 'off');
% varargout = sigstar({[1 2]}, [early_vs_late_pval]);
% 
% title_str = ['early: ' num2str(early_mean) ', late: ' num2str(late_mean) ...
%     ' ranksum: ' num2str(early_vs_late_pval) ' num targ: ' num2str(num_consistent_targets)];
% title(title_str);
% 
% %%
% export_fig(h, fullfile(save_dir, [data_str ' early_late_over_T' '.eps']));
% 
% 
% end