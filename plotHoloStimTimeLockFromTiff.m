function plotHoloStimTimeLockFromTiff(fname, mask, voltageRec, wd)

%{
Function to plot the response of neurons to the holo_stim
fname => tiff file
mask => mask of neurons, simple N*M matrix obtained with bwlabel
voltageRec => data from the voltage recording
wd => window to plot before/after neuron
%}
    if nargin < 4
        wd = 100;
    end

    % to obtain the timing of the holostims
    holoStim = voltageRec(:,8);
    meanHoloStim = nanmean(holoStim);
    holoStim(holoStim<meanHoloStim) = 0;
    holoStim(holoStim>meanHoloStim) = 1;
    [~,locs] = findpeaks(holoStim);

    % to obtain number of subplots
    numNeuron = max(max(mask));
    numSubp = ceil(sqrt(numNeuron));
    figure();

    % to obtain the activity of the neuron
    info = imfinfo(fname);
    num_images = numel(info);
    for k = 1:num_images
        A = imread(fname, k);
        for nn = 1:numNeuron
            auxIndex = find(locs(nn)<=temp, 1, 'first');
            minIndex = max(auxIndex - wd, 1);
            maxIndex = min(auxIndex + wd, size(botData,1));
            subplot(numSubp,numSubp,nn-1)
            plot(temp(minIndex:maxIndex), botData(minIndex:maxIndex, nn))
            vline(temp(auxIndex), 'r-')
            title(['neuron: ', int2str(nn)])
        end

    end
end