function plotNeuronsBaseline(baseActivity, CComp, YrA, totalNeurons)
%{
Function to plot the temporal activity of neurons collected during Baseline
to select the best neurons.
baseActivity -> activity during baseline
CComp -> C_on from holostim period given by onacid
YrA -> background noise of C
totalNeurons -> amount of neurons to be displayed

%}
    if nargin < 4
        totalNeurons = 20;
    end
	Sm = nanstd(baseActivity(1:end,10:end),2)./nanmean(baseActivity(1:end,10:end),2);
    S = nanstd(baseActivity(1:end,10:end),2);
	[~, indm] = sort(Sm, 'descend');
    [~, ind] = sort(S, 'descend');
    disp('Neurons from best to worst Sm: \n');
	indm(1:totalNeurons) 
    disp('Neurons from best to worst S: \n');
	ind(1:totalNeurons)
    % plot std/mean
	figure()
    sgtitle('Std/mean')
    for idx=1:totalNeurons
		subplot(4,5,idx)
		plot(baseActivity(ind(idx), :)');
		title(['ROI ' int2str(ind(idx))]);
    end
    % plot std
    figure()
    sgtitle('Std')
    for idx=1:totalNeurons
		subplot(4,5,idx)
		plot(baseActivity(ind(idx), :)');
		title(['ROI ' int2str(ind(idx))]);
    end
    % plot C and Cnoise
    CNoise = CComp + YrA;
    figure()
    sgtitle('Std')
    for idx=1:totalNeurons
		subplot(4,5,idx)
		plot(CNoise(ind(idx), :)');
        hold on
        plot(CComp(ind(idx), :)');
		title(['ROI ' int2str(ind(idx))]);
    end
end
    