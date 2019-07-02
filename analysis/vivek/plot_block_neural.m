function plot_block_neural(block_str, block_data, event_struct, frame_rate, win_zoom, save_dir, save_bool, plot_big)
%Plots neural activity progression in a block (base, pretrain, bmi)
%block_data: neural data from 'base_bmiAct2cursor.m'
%event_struct: 2d struct array, num_blocks x num_events t
%   fields: data, label, valid
%frame_rate: frames / sec
%win_zoom: [start stop] in units of minutes

time_scale  =1/(frame_rate*60); %frames per minute
%%
%Plot each block's activity.  Plot f and f0 overlaid, plot dff
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 

%F, F0
%dff
%dff_z
%For each: 
%   event, zoom event

for i=1:length(block_data)
    E_id = block_data(i).bData.E_id; 
    num_BMI = length(E_id);    

%     %1) F F0
%     f   = block_data(i).data_valid.bmiAct; 
%     f0  = block_data(i).data_valid.baseVector;
%     f_f0_name = [block_str num2str(i) '_f_f0'];
%     t   = (1:length(f))*time_scale; 
%     
%     label = f_f0_name; 
%     [h_cell] = ...
%         plot_F_F0(E_id, t, f, f0, label, event_struct(i,:), save_dir, save_bool, E_color, plot_big);
%     
%     %1.1) F F0 Zoom
%     zoom_samples = find(t>=win_zoom(1) & t<=win_zoom(2));
%     t_zoom = t(zoom_samples); 
%     f_zoom = f(:,zoom_samples); 
%     f0_zoom = f0(:,zoom_samples);     
%     zoom_event_struct = zoom_event(event_struct(i,:), win_zoom); 
%     
%     label = [f_f0_name '_zoom']; 
%     [h_cell] = ...
%         plot_F_F0(E_id, t_zoom, f_zoom, f0_zoom, label, zoom_event_struct, save_dir, save_bool, E_color, plot_big);   
% 
%     %2) dff
%     dff_name    = [block_str num2str(i) '_dff']; 
%     dff         = block_data(i).dff.'; 
%     t           = (1:length(dff))*time_scale;
%     
%     label       = dff_name; 
%     [h_cell] = plot_n_0(E_id, t, dff, label, event_struct, save_dir, save_bool, E_color, plot_big);
%     
%     %2.1) dff zoom
%     zoom_samples = find(t>=win_zoom(1) & t<=win_zoom(2));
%     t_zoom = t(zoom_samples); 
%     dff_zoom = dff(zoom_samples,:);   
%     zoom_event_struct = zoom_event(event_struct(i,:), win_zoom); 
%     
%     label = [dff_name '_zoom']; 
%     [h_cell] = plot_n_0(E_id, t_zoom, dff_zoom, label, zoom_event_struct, save_dir, save_bool, E_color, plot_big);       
% 
%     %3) dffz
    dffz_name    = [block_str num2str(i) '_dffz']; 
    dff_z         = block_data(i).dff_z.'; 
    t           = (1:length(dff_z))*time_scale;
    
    label       = dffz_name; 
    [h_cell] = plot_n_0(E_id, t, dff_z, label, event_struct(i,:), save_dir, save_bool, E_color, plot_big);
    
    %3.1) dffz zoom
    zoom_samples = find(t>=win_zoom(1) & t<=win_zoom(2));
    t_zoom = t(zoom_samples); 
    dffz_zoom = dff_z(zoom_samples,:);   
    zoom_event_struct = zoom_event(event_struct(i,:), win_zoom); 
    
    label = [dffz_name '_zoom']; 
    [h_cell] = plot_n_0(E_id, t_zoom, dffz_zoom, label, zoom_event_struct, save_dir, save_bool, E_color, plot_big);  
end

end

function zoom_event_struct = zoom_event(event_struct, zoom_win)
    zoom_event_struct = event_struct; 
    for i = 1:length(event_struct)
        sel_idxs = ...
            find((event_struct(i).data >= zoom_win(1)) & ...
            (event_struct(i).data <= zoom_win(2))); 
        zoom_event_struct(i).data = zoom_event_struct(i).data(sel_idxs);
    end
end

function [h_cell] = plot_F_F0(E_id, t, f, f0, label, event_struct, save_dir, save_bool, E_color, plot_big)
%Inputs: 
%num_BMI - number of channels to plot
% t - 1 x num_timesamples
%f - fluorescence
%f0 - baseline fluorescence
%label - name to save plot
%event_struct - struct with events to plot with vline
%   fields: data, label, valid
    %
    h_cell = {}; 
    %1) F F0
    plot_name = label; 
    num_BMI = length(E_id); 
    n_cell = cell(num_BMI, 1); 
    t_cell = cell(num_BMI, 1); 
    for j =1:num_BMI
        n_cell{j} = [f(j,:)' f0(j,:)']; 
        t_cell{j} = t; 
    end
    [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color);
    set(gca,'TickDir','out');
    xlabel('time (min)'); 
    title(label); 
    if(plot_big)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    end
    if(save_bool)
        export_fig(h, fullfile(save_dir, [plot_name '.eps'])); 
        export_fig(h, fullfile(save_dir, [plot_name '.png']));
    end    
    h_cell{1} = h; 
    
    
    %2) loop events: 
    num_event_types = length(event_struct); 
    if(num_event_types > 1)
        for event_i = 1:num_event_types
            if(event_struct(event_i).valid)
                plot_name = [label '_' event_struct(event_i).label]; 
                [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color);
                vline(event_struct(event_i).data); 
                set(gca,'TickDir','out');
                xlabel('time (min)'); 
                title(event_struct(event_i).label);
                if(plot_big)
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                end
                if(save_bool)
                    export_fig(h, fullfile(save_dir, [plot_name '.eps'])); 
                    export_fig(h, fullfile(save_dir, [plot_name '.png']));
                end
                h_cell{end+1} = h; 
            end        
        end
    end
        
    %Plot all events together: 
    plot_name = [label '_all_events']; 
    [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color);    
    num_event_types = length(event_struct); 

    v_h = []
    event_label_cell = {}
    for event_i = 1:num_event_types
        if(event_struct(event_i).valid)
            h_i = vline(event_struct(event_i).data, event_struct(event_i).line); 
            v_h = [v_h h_i(1)]; 
            event_label_cell{event_i} = event_struct(event_i).label; 
        end        
    end
    legend(v_h, event_label_cell); 
    h_cell{end+1} = h;
    set(gca,'TickDir','out');
    xlabel('time (min)'); 
    title('all events');
    if(plot_big)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    end
    if(save_bool)
        export_fig(h, fullfile(save_dir, [plot_name '.eps'])); 
        export_fig(h, fullfile(save_dir, [plot_name '.png']));
    end
    
end

function [h_cell] = plot_n_0(E_id, t, n, label, event_struct, save_dir, save_bool, E_color, plot_big)
%Plots time series 'n' and overlays a plot of horizontal line 0
%     plot_big = 0; 

    h_cell = {}; 
    num_BMI = length(E_id);

    n_cell = cell(num_BMI, 1); 
    t_cell = cell(num_BMI, 1); 
    for j =1:num_BMI
        y_j = n(:,j);
        n_cell{j} = [y_j zeros(size(y_j))]; 
        t_cell{j} = t;
    end
    
    %1) dff 0
    plot_name = label; 
    [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color);
    set(gca,'TickDir','out');
    xlabel('time (min)'); 
    title(label); 
    if(plot_big)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    end
    if(save_bool)
        export_fig(h, fullfile(save_dir, [plot_name '.eps'])); 
        export_fig(h, fullfile(save_dir, [plot_name '.png']));
    end    
    h_cell{end+1} = h;  
    
    %2) loop events: 
    num_event_types = length(event_struct); 
    if(num_event_types > 1)
        for event_i = 1:num_event_types
            
            if(event_struct(event_i).valid)
                plot_name = [label '_' event_struct(event_i).label]; 
                [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color);
                vline(event_struct(event_i).data); 
                set(gca,'TickDir','out');
                xlabel('time (min)'); 
                title([label event_struct(event_i).label]);
                if(plot_big)
                    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
                end
                if(save_bool)
                    export_fig(h, fullfile(save_dir, [plot_name '.eps'])); 
                    export_fig(h, fullfile(save_dir, [plot_name '.png']));
                end
                h_cell{end+1} = h; 
            end
            
        end
    end
    
    %Plot all events together: 
    plot_name = [label ' all events']; 
    [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color);    
    num_event_types = length(event_struct); 

    v_h = []
    event_label_cell = {}
    for event_i = 1:num_event_types
        if(event_struct(event_i).valid)
            h_i = vline(event_struct(event_i).data, event_struct(event_i).line); 
            v_h = [v_h h_i(1)]; 
            event_label_cell{event_i} = event_struct(event_i).label; 
        end        
    end
    legend(v_h, event_label_cell); 
    h_cell{end+1} = h;
    set(gca,'TickDir','out');
    xlabel('time (min)'); 
    title([label ' all events']);
    if(plot_big)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    end
    if(save_bool)
        export_fig(h, fullfile(save_dir, [plot_name '.eps'])); 
        export_fig(h, fullfile(save_dir, [plot_name '.png']));
    end    
    
end
    

%     %2) F F0 zoom
%     zoom_samples = find(t>=win_zoom(1) & t<=win_zoom(2));
%     t_zoom = t(zoom_samples); 
%     f_zoom = f(:,zoom_samples); 
%     f0_zoom = f0(:,zoom_samples); 
%     
%     
%     %3) 
%     
%     
%     
%     
%     
%     %4) dff
%     plot_name = ['pretrain' num2str(i) '_dff_stim']; 
%     n = block_data(i).dff.';
%     n_cell = cell(num_BMI, 1); 
%     t_cell = cell(num_BMI, 1); 
%     for j =1:num_BMI
%         y_j = n(:,j);
%         n_cell{j} = [y_j zeros(size(y_j))]; 
%         t_cell{j} = (1:length(y_j))*time_scale;
%     end    
%     [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color)
%     title('pretrain dff'); 
%     xlabel('time (min)');     
%     vline(stim_time); 
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     export_fig(h, fullfile(plot_dir, [plot_name '.eps'])); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.png']));       
%     
%     %5) dff zoom
%     time_zoom = 5; 
%     num_samples = 5*30*60;
%     stim_zoom = stim_time(stim_time<=time_zoom);    
%     
%     plot_name = ['pretrain' num2str(i) '_dff_stim_zoom' num2str(time_zoom)]; 
%     n = block_data(i).dff.';
%     n_cell = cell(num_BMI, 1); 
%     t_cell = cell(num_BMI, 1); 
%     for j =1:num_BMI
%         y_j = n(1:num_samples,j);
%         n_cell{j} = [y_j zeros(size(y_j))]; 
%         t_cell{j} = (1:length(y_j))*time_scale;
%     end    
%     [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color)
%     title('pretrain dff'); 
%     xlabel('time (min)');     
%     vline(stim_zoom); 
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     export_fig(h, fullfile(plot_dir, [plot_name '.eps'])); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.png']));  
%     
%     %6) dff zscore
%     plot_name = ['pretrain' num2str(i) '_dffz_stim']; 
%     n = zscore(block_data(i).dff.');
%     n_cell = cell(num_BMI, 1); 
%     t_cell = cell(num_BMI, 1); 
%     for j =1:num_BMI
%         y_j = n(:,j);
%         n_cell{j} = [y_j zeros(size(y_j))]; 
%         t_cell{j} = (1:length(y_j))*time_scale;
%     end    
%     [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color)
%     title('pretrain dff zscore'); 
%     xlabel('time (min)');     
%     vline(stim_time); 
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     export_fig(h, fullfile(plot_dir, [plot_name '.eps'])); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.png']));     
%     
%     %7) dff zscore zoom    
%     time_zoom = 5; 
%     num_samples = 5*30*60;
%     stim_zoom = stim_time(stim_time<=time_zoom);     
%     
%     plot_name = ['pretrain' num2str(i) '_dffz_stim_zoom' num2str(time_zoom)]; 
%     n = zscore(block_data(i).dff.');
%     n_cell = cell(num_BMI, 1); 
%     t_cell = cell(num_BMI, 1); 
%     for j =1:num_BMI
%         y_j = n(1:num_samples,j);
%         n_cell{j} = [y_j zeros(size(y_j))]; 
%         t_cell{j} = (1:length(y_j))*time_scale;
%     end    
%     [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color)
%     title('pretrain dff zscore zoom'); 
%     xlabel('time (min)');     
%     vline(stim_zoom); 
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     export_fig(h, fullfile(plot_dir, [plot_name '.eps'])); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.png']));        
% end
% 
% end




    
    
%     %2) F F0 Stim
%     valid_idxs = block_data(i).data_valid.valid_idxs;
%     event_time = time_scale*find(pretrain(i).data.holoDelivery(valid_idxs)); 
% 
%     plot_name = ['pretrain' num2str(i) '_f_f0_stim']; 
% %     [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color)
%     vline(stim_time); 
%     xlabel('time (min)'); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.eps'])); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.png']));  
% 
%     %3) Zoom F F0 Stim
%     time_zoom = 5; 
%     num_samples = 5*30*60;
%     stim_zoom = stim_time(stim_time<=time_zoom);
% 
%     plot_name = ['pretrain' num2str(i) '_f_f0_stim_zoom' num2str(time_zoom)]; 
%     f   = block_data.data_valid.bmiAct; 
%     f0  = block_data.data_valid.baseVector;
%     n_cell = cell(num_BMI, 1); 
%     t_cell = cell(num_BMI, 1); 
%     for i =1:num_BMI
%         n_cell{i} = [f(i,1:num_samples)' f0(i,1:num_samples)']; 
%         t_cell{i} = (1:length(f(i,1:num_samples)))*time_scale;
%     end
%     [h, offset_vec] = plot_E_activity_mult_trace(t_cell, n_cell, E_id, E_color)
%     vline(stim_zoom); 
%     xlabel('time (min)'); 
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     export_fig(h, fullfile(plot_dir, [plot_name '.eps'])); 
%     export_fig(h, fullfile(plot_dir, [plot_name '.png'])); 
% 
% end