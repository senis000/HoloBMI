function plotNeuronsEnsemble(baseActivity, ensembleNeurons)
%{
Function to plot the temporal activity of neurons collected during Baseline
to select the best neurons.
baseActivity -> activity during baseline
CComp -> C_on from holostim period given by onacid
YrA -> background noise of C
totalNeurons -> amount of neurons to be displayed

%}
    totalNeurons = length(ensembleNeurons);
    
    subplotnmb = ceil(totalNeurons/2);
	figure('Position', [300,300, subplotnmb*200, 400])
    %sgtitle('Std/mean')
    for idx = 1:length(ensembleNeurons)
		subplot(2,subplotnmb,idx)
		plot(baseActivity(ensembleNeurons(idx), :)');
		title(['ROI ' int2str(ensembleNeurons(idx))]);
    end
end
    