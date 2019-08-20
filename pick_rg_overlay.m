function [rg_struct, num_overlays] = pick_rg_overlay(red_im, green_im, rg_struct, num_overlays)
%Output: 
%rg_struct fields: 
% rg_struct(num_overlays).im = rg; 
% rg_struct(num_overlays).rg_minmax_perc = param_vec; 
% rg_struct(num_overlays).rg_minmax = [r_min r_max g_min g_max]; 
% rg_struct(num_overlays).r_min = r_min; 
% rg_struct(num_overlays).r_max = r_max; 
% rg_struct(num_overlays).r_min_perc = r_min_perc; 
% rg_struct(num_overlays).r_max_perc = r_max_perc; 
% rg_struct(num_overlays).g_min = g_min; 
% rg_struct(num_overlays).g_max = r_max; 
% rg_struct(num_overlays).g_min_perc = g_min_perc; 
% rg_struct(num_overlays).g_max_perc = g_max_perc; 

screen_size = get(0,'ScreenSize');
%%
%%Prompt user to change g_min_perc, g_max_perc, r_min_perc, r_max_perc:
overlay_complete_bool = 0;
param_vec = [0 100 0 100]; 
%initialize parameters: 
while(~overlay_complete_bool)
    disp('Current [r_min_perc r_max_perc g_min_perc g_max_perc ]:')
    disp(param_vec); 
    param_vec = input('enter [r_min_perc r_max_perc g_min_perc g_max_perc]: ');
    while(length(param_vec) ~= 4)
        disp('Error in input!'); 
        param_vec = input('enter [ r_min_perc r_max_perc g_min_perc g_max_perc]: ');
    end
%     [g_min_perc, g_max_perc, r_min_perc, r_max_perc] = param_vec; 
    r_min_perc = param_vec(1); 
    r_max_perc = param_vec(2); 
    g_min_perc = param_vec(3); 
    g_max_perc = param_vec(4); 

    
    [green_s, g_min, g_max] = scale_im(green_im, g_min_perc, g_max_perc);
    [red_s, r_min, r_max]   = scale_im(red_im, r_min_perc, r_max_perc);
    
    %----------------------------------------------------------------------
    rg = zeros(size(green_im,1), size(green_im,2), 3); 
    rg(:,:,1) = red_s;
    rg(:,:,2) = green_s; 
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(rg); 
    axis square
    title(['rmin: ' num2str(r_min_perc) ' rmax: ' num2str(r_max_perc) ' gmin: ' num2str(g_min_perc) ' gmax: ' num2str(g_max_perc)]);     
    
    %----------------------------------------------------------------------
    in = input('Want to Keep Overlay? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'y'))    
        %Check you don't already have it: 
        already_added = 0; 
        for i = 1:num_overlays
            if(sum(rg_struct(i).rg_minmax ~= param_vec) ==0)
                already_added = 1; 
                disp('already saved it!'); 
                break;
            end
        end
        
        if(~already_added)
            num_overlays = num_overlays+1; 
            rg_struct(num_overlays).im = rg; 
            rg_struct(num_overlays).rg_minmax_perc = param_vec; 
            rg_struct(num_overlays).rg_minmax = [r_min r_max g_min g_max]; 
            rg_struct(num_overlays).r_min = r_min; 
            rg_struct(num_overlays).r_max = r_max; 
            rg_struct(num_overlays).r_min_perc = r_min_perc; 
            rg_struct(num_overlays).r_max_perc = r_max_perc; 

            rg_struct(num_overlays).g_min = g_min; 
            rg_struct(num_overlays).g_max = r_max; 
            rg_struct(num_overlays).g_min_perc = g_min_perc; 
            rg_struct(num_overlays).g_max_perc = g_max_perc;         
        end
    end
    
    in = input('Want to Add Overlays? y/n:   ', 's');
    in = lower(in);    
    if(strcmp(in, 'n'))
        overlay_complete_bool = 1;
    end
end