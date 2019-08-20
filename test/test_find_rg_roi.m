
%%
g_path = fullfile('E:\ines_e\redgreen-000', 'AVG_green.tif'); 
r_path = fullfile('E:\ines_e\redgreen-000', 'AVG_red.tif');

%%
g_mean = imread(g_path); 
r_mean = imread(r_path);   

%%
h = figure;
imagesc(r_mean);
colormap bone
axis square
title('r mean'); 

%
h = figure;
imagesc(g_mean);
colormap bone
axis square
title('g mean'); 

%%
%Initialize data structure to save overlays:
rg_struct = struct(...
    'im', [], ...
    'rg_minmax_perc', [], ...
    'r_min', [], ...
    'r_min_perc', [], ...
    'g_min', [], ...
    'g_min_perc', [], ...
    'r_max', [], ...
    'r_max_perc', [], ...
    'g_max', [], ...
    'g_max_perc', []...
    ); 
i = 0; 
%%
g_min_perc = 0.1;
g_max_perc = 100.;%99.99; 
% min_perc = 0.01;
% max_perc = 99.99; 
[green_s, g_min, g_max] = scale_im(g_mean, g_min_perc, g_max_perc);
% h = figure;
% imagesc(green_s); colormap bone 
% axis square
%
r_min_perc = 0;
r_max_perc = 99.9%%99.5; %99.99; 

[red_s, r_min, r_max] = scale_im(r_mean, r_min_perc, r_max_perc);
% h = figure;
% imagesc(red_s); colormap bone 
% axis square
% function im_s = scale_im(im, min_perc, max_perc)

rg = zeros(size(g_mean,1), size(g_mean,2), 3); 
rg(:,:,1) = red_s;
rg(:,:,2) = green_s; 
h = figure;
imagesc(rg); 
axis square
title(['rmin: ' num2str(r_min_perc) ' rmax: ' num2str(r_max_perc) ' gmin: ' num2str(g_min_perc) ' gmax: ' num2str(g_max_perc)]); 

%%
%if we want the image saved: 
i = i+1;
rg_struct(i).im = rg; 

rg_struct(i).r_min = r_min; 
rg_struct(i).r_min_perc = r_min_perc; 
rg_struct(i).g_min = g_min; 
rg_struct(i).g_min_perc = g_min_perc; 

rg_struct(i).r_max = r_max; 
rg_struct(i).r_max_perc = r_max_perc; 
rg_struct(i).g_max = g_max; 
rg_struct(i).g_max_perc = g_max_perc; 

%%
r = input('enter vector')

%%
%The image we want to select ROI's on.  
rg_sel = rg_struct(i)
%Liked: r_max: 99.9

%%
h = figure;
imagesc(rg_sel.im); axis square;

%%
mask = zeros(512,512); 
mat = rg_sel.im; 
temp = mask;
dim=8;
[mask, x_cen, y_cent] = addcell (mat,temp,dim);

%%
h = figure;
imagesc(mask)

%%
h = figure;
imagesc(mat); 
axis square
hold on;
scatter(x_cen, y_cent, 'filled', 'r'); 

%%
toplot = 1
[x,y]=findCenter(mask, mat, toplot)


%%
Im = mat;
[x,y] = findCenter (mask, Im, false);
h = figure;
imagesc(Im), colormap bone, caxis([-0 nanmean(nanmean(Im(:)))*4]), hold on, scatter (x,fliplr(y), 'filled', 'r'), hold off, axis square

%%
numArea = 1;
[rois] = test_deleteMask(Im, mask, numArea);

%%
%Implement for choosing ROI, then re-use for adding two channels of ROI: 
%Show the image, with roi's scattered on it.  Press 'y' to add roi, 'n' to
%finish.  
%Draw ROI
%Show background image with just the new ROI on it
%press 'y' to keep the ROI, press 'n' to delete it
%

%%
% Show the red and green channels:
close all;
%GREEN
im_plot = g_mean; 
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(im_plot); axis square; colormap bone; title('Green Mean'); 
%RED
im_plot = r_mean;
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(im_plot); axis square; colormap bone; title('Red Mean'); 

%
disp('Adding ROIs to image!'); 
roi_init_bool = 1; 
if(roi_init_bool)
    Im_roi = Im; 
    num_row = size(Im_roi,1); 
    num_col = size(Im_roi,2); 
    roi_mask        = zeros(num_row, num_col); 
    roi_mask_bin    = zeros(num_row, num_col); 
    roi_bin_cell    = {}; 
    num_rois = 0; 
end

roi_complete_bool = 0; 
screen_size = get(0,'ScreenSize');
while(~roi_complete_bool)
    disp('Current Roi Image: '); 
    if(sum(ismember(findall(0,'type','figure'),h0)))
        close(h0)
    end
    if(sum(ismember(findall(0,'type','figure'),h1)))
        close(h1)
    end
    h0 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(Im_roi); caxis([-0 nanmean(nanmean(Im(:)))*4]); axis square
    title(['Num ROIs added: ' num2str(num_rois) '  Add ROI? y/n']); 
    in = input('Want to Add ROI? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'y'))
        disp('Draw the ROI (click, hold click, draw on image)...')
        title('Draw ROI'); 
        drawRois(numArea);  %draw rois

        %Collect roi data: 
        rois = [];
        hF = gcf;
        hP = findobj(hF, 'Tag', 'ROIPatch'); %finds rois
        for rr = 1:size(hP, 1)
            rois(:,:,rr) = getUD(hP(rr), 'binroi'); %gets rois
        end

        %Add single ROI 
        Im_roi_i = Im;
        Im_roi_i(:,:,3) = rois; 
        h1 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
        imagesc(Im_roi_i); axis square; title(['Added Candidate Roi # ' num2str(num_rois+1) ' Keep? y/n']);
        %Check whether to keep ROI:
        in = input('Keep ROI? y/n:   ', 's');
        in = lower(in);
        if(isempty(in) || strcmp(in, 'y'))
            num_rois = num_rois+1; 
            %This ROI mask binary:
            roi_bin_cell{num_rois}      = rois;

            roi_mask(find(rois))        = num_rois;
            roi_mask_bin(find(rois))    = 1;   
            Im_roi(:,:,3)               = roi_mask_bin;
        end
    else
        roi_complete_bool = 1;
        disp('done')
        if(sum(ismember(findall(0,'type','figure'),h0)))
            close(h0)
        end
        if(sum(ismember(findall(0,'type','figure'),h1)))
            close(h1)
        end
        h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
        imagesc(Im_roi); caxis([-0 nanmean(nanmean(Im(:)))*4]); axis square
        title(['ROI addition complete! Num ROIs: ' num2str(num_rois)]);         
    end
end

%%
h = figure;
imagesc(Im); axis square

%%
h = figure;
im_clim = [-0 nanmean(nanmean(Im(:)))*1000];
imagesc(Im); axis square; caxis(im_clim);

%%
%test function decomposition: 
plot_images = struct('im', [], 'label', ''); 
plot_images(1).im = g_mean; 
plot_images(1).label = 'Green Mean'; 
plot_images(2).im = r_mean; 
plot_images(2).label = 'Red Mean'; 

%Initialize roi_data: 
roi_data.im_bg              = Im; %needs to be rgb, {num_row x num_col x 3}
roi_data.im_roi             = Im;
[roi_data.num_row, roi_data.num_col, roi_data.num_color]    = size(Im_roi); 
roi_data.roi_mask            = zeros(roi_data.num_row, roi_data.num_col); 
roi_data.roi_mask_bin        = zeros(roi_data.num_row, roi_data.num_col); 
roi_data.roi_bin_cell        = {}; 
roi_data.num_rois           = 0; 

%drawcell

%%

h = figure;
imagesc(roi_data.roi_mask_bin); 
%%
[roi_data] = drawcell(plot_images, roi_data);

%%
h = figure;
imagesc(roi_data.roi_mask); axis square; colorbar;

%%
h = figure;
imagesc(roi_mask); 

%%
% roi_data.im_roi     = Im; 
% 
% 
% 
% roi_init_bool = 1; 
% if(roi_init_bool)
%     Im_roi = Im; 
%     [num_row = size(Im_roi,1); 
%     num_col = size(Im_roi,2); 
%     roi_mask        = zeros(num_row, num_col); 
%     roi_mask_bin    = zeros(num_row, num_col); 
%     roi_bin_cell    = {}; 
%     num_rois = 0; 
% end



%%
%Test delete ROI works...
%(Can improve delete ROI)

%%
% select areas to delete

numArea = 1;
drawRois(numArea);  %draw rois
for ind = 1:(numArea+3)
    waitforbuttonpress;  %pauses until finished
end
rois = [];
hF = gcf;
hP = findobj(hF, 'Tag', 'ROIPatch'); %finds rois
for rr = 1:size(hP, 1)
    rois(:,:,rr) = getUD(hP(rr), 'binroi'); %gets rois
end
todel = unique(nansum(rois,3).*mask); %finds index of rois