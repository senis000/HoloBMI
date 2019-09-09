function [roi_data] = draw_roi(plot_images, roi_data)
%Allows user to draw shapes onto an image

%%
screen_size = get(0,'ScreenSize');
% Show the red and green channels:
close all;
for plot_i=1:length(plot_images)
    im_plot = plot_images(plot_i).im; 
    im_title = plot_images(plot_i).label; 
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(im_plot); axis square; colormap bone; title(im_title);    
end

disp('Adding ROIs to image!'); 
roi_complete_bool = 0; 
while(~roi_complete_bool)
    %Close images:
    if(exist('h0'))
        if(sum(ismember(findall(0,'type','figure'),h0)))
            close(h0)
        end
    end
    if(exist('h1'))
        if(sum(ismember(findall(0,'type','figure'),h1)))
            close(h1)
        end
    end
    
    disp(['Current Roi Image, Num Rois: ' num2str(roi_data.num_rois)]); 
    h0 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    imagesc(roi_data.im_roi); axis square
    title(['Num ROIs added: ' num2str(roi_data.num_rois) '  Add ROI? y/n']); 
    in = input('Want to Add ROI? y/n:   ', 's');
    in = lower(in);
    if(isempty(in) || strcmp(in, 'y'))
        disp('Draw the ROI (click, hold click, draw on image)...')
        title('Draw ROI'); 
        numArea = 1; %Only draw one ROI
        drawRois(numArea);  %draw rois

        %Collect roi data: 
        rois = [];
        hF = gcf;
        hP = findobj(hF, 'Tag', 'ROIPatch'); %finds rois
        for rr = 1:size(hP, 1)
            rois(:,:,rr) = getUD(hP(rr), 'binroi'); %gets rois
        end

        %Add single ROI 
        Im_roi_i = roi_data.im_bg;
        Im_roi_i(:,:,3) = rois; 
        h1 = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
        imagesc(Im_roi_i); axis square; title(['Added Candidate Roi # ' num2str(roi_data.num_rois+1) ' Keep? y/n']);
        %Check whether to keep ROI:
        in = input('Keep ROI? y/n:   ', 's');
        in = lower(in);
        if(isempty(in) || strcmp(in, 'y'))
            roi_data.num_rois = roi_data.num_rois+1; 
            %This ROI mask binary:
            roi_data.roi_bin_cell{roi_data.num_rois}      = rois;

            roi_data.roi_mask(find(rois))        = roi_data.num_rois;
            roi_data.roi_mask_bin(find(rois))    = 1;   
            roi_data.im_roi(:,:,3)               = roi_data.roi_mask_bin;
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
        imagesc(roi_data.im_roi); axis square
        title(['ROI addition complete! Num ROIs: ' num2str(roi_data.num_rois)]);         
    end
end