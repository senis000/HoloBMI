function [rois] = test_deleteMask(Im, mask, numArea)
%{ 
Function to delete neurons of a mask 
Im -> Raw image
Mask -> mask from where we want to delete neurons
numArea -> amount of areas from where to delete

An image will open, you can maximize the image and then draw rois around
the neurons that wish to be deleted. When finished press "enter"
be careful not to include neurons you don't want to delete
be careful not to overlap rois
you can select some more areas that the number defined in numArea, but not
much more. If you select less you will have to press more times "enter" at
the end of the selection of areas
%}

    if nargin < 3
        numArea = 1;
    end
    close all
    % display the image and mask
    [x,y] = findCenter (mask, Im, false);
    imagesc(Im), colormap bone, caxis([-0 nanmean(nanmean(Im(:)))*4]), hold on, scatter (x,fliplr(y), 'filled', 'r'), hold off, axis square
    
    % select areas to delete
    
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
    
%     todel = unique(nansum(rois,3).*mask); %finds index of rois
%     for tt=1:length(todel)
%         mask(mask==todel(tt)) = 0;  %delete neurons 
%     end
%     
%     % reorder the mask
%     mask(mask>0) = 1;
%     mask = bwlabel(mask);
%     
%     %plot for sanity check
%     findCenter (mask, Im);  
end

