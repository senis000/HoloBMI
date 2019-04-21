function [data_sem, data_mean, rho, pval, h] = ...
    plot_data_and_corr_over_cond(num_var, num_cond, data_over_cond_var, data_valid, x_data, data_str, save_dir, ...
    x_win_bool, x_win, y_win_bool, y_win)
%INPUT:
%data_over_cond_var:
% data matrix:
% size: num_cond x num_variables

    MarkerSize  = 1; %20
%Just Data + Error Bar:
    [data_sem, data_mean] = sem(data_over_cond_var.', data_valid.');
    h = figure;
    hold on;
%    hErr = errorbar(x_data, data_mean, data_sem, '.-', 'MarkerSize', MarkerSize, 'color', [0 0 0], 'LineWidth', 2);    
     errbar(x_data, data_mean, data_sem, 'LineWidth', 2, 'color', [0 0 0]);
     plot(x_data, data_mean, '.-', 'MarkerSize', 20, 'color', [0 0 0], 'LineWidth', 2);
    
    
    %x_data = (1:num_var)';
    x_data  = x_data(:);
    y_data  = data_mean(:);
    
    invalid_idxs     = find(isnan(y_data));
    
    x_data(invalid_idxs) = [];
    y_data(invalid_idxs) = [];
    
    [rho, pval] = corr(x_data, y_data);
    title(['rho: ' num2str(rho) ' pval: ' num2str(pval)]);
    
    xmin    = x_data(1)-1;
    xmax    = x_data(end)+1;
    
    if(~x_win_bool)
        xlim([xmin xmax]);
    else
        xlim(x_win);
    end
    
    if(y_win_bool)
        ylim(y_win);
    else
%     ylim([25 85]); 
%     ylim([-1 1.5]);
%       ylim([0 5]);
%       ylim([-1 2]);        
    end
    


%align:
% ylim([-1.5 1.5]);

    if(num_var > 7)
        if(ceil(num_var/2) > 1 && ceil(num_var/2) < num_var)
            set(gca, 'XTick', [1 ceil(num_var/2) num_var]);  
        else
            set(gca, 'XTick', [1 num_var]);  
        end
    else
        set(gca, 'XTick', [1:num_var]);  
    end
    
%     ymin = min(data_mean-data_sem);
%     ymax = max(data_mean+data_sem);
%     set(gca, 'YTick', [ymin (ymax-ymin)/2 ymax]);
    
    %set(gca, 'XTick', 1:2:max_num_epochs);
    
%    ylim([0.6 3]);

    
    set(gca, 'TickDir', 'out');
    set(gca, 'box', 'off');
    
    export_fig(h, fullfile(save_dir, [data_str '_errorbar.eps']));