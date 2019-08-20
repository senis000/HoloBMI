function plotSelNeuronsBaseline(baseActivity, CComp, YrA, totalNeurons, sel)
%{
Function to plot the temporal activity of neurons collected during Baseline
to select the best neurons.
baseActivity -> activity during baseline
CComp -> C_on from holostim period given by onacid
YrA -> background noise of C
totalNeurons -> amount of neurons to be displayed
sel -> idxs of the baseActivity to plot

%}
    if nargin < 4
        totalNeurons = 20;
    end
    totalNeurons = min(totalNeurons, length(sel)); 
    
    subplotnmb = ceil(sqrt(totalNeurons));
    
    prefix_discard = 10; 
	Sm = nanstd(baseActivity(sel,prefix_discard:end),0,2)./nanmean(baseActivity(sel,prefix_discard:end),2);
    S = nanstd(baseActivity(sel,prefix_discard:end),0,2);
	[~, indm] = sort(Sm, 'descend');
    [~, ind] = sort(S, 'descend');
    disp('Neurons from best to worst Sm: \n');
	sel(indm(1:totalNeurons))
    disp('Neurons from best to worst S: \n');
	sel(ind(1:totalNeurons))
    % plot std/mean
	figure()
    %sgtitle('Std/mean')
    for idx=1:totalNeurons
		subplot(subplotnmb,subplotnmb,idx)
		plot(baseActivity(sel(ind(idx)), :)');
		title(['ROI ' int2str(sel(ind(idx)))]);
    end
    % plot std
%     figure()
% %    sgtitle('Std')
%     for idx=1:totalNeurons
% 		subplot(4,5,idx)
% 		plot(baseActivity(indm(idx), :)');
% 		title(['ROI ' int2str(indm(idx))]);
%     end

%     % plot C and Cnoise
%     if(~isempty(CComp))
%         CNoise = CComp + YrA;
%         figure()
%     %    sgtitle('HoloStim')
%         for idx=1:totalNeurons
%             subplot(subplotnmb,subplotnmb,idx)
%             plot(CNoise(indm(idx), :)');
%             hold on
%             plot(CComp(indm(idx), :)');
%             title(['ROI ' int2str(indm(idx))]);
%         end
%     end
end
    