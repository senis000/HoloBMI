function plotHoloStimTimeLockFromTiff(fname, mask, voltageRec, wd)

%{
Function to plot the response of neurons to the holo_stim
fname => tiff file
mask => mask of neurons, simple N*M matrix obtained with bwlabel
voltageRec => data from the voltage recording
wd => window to plot before/after neuron

TODO CHECK FRAMES!!!! WHICH VOLTAGE OUTPUT IS IT
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

    % to load the image info one by one 
    info = imfinfo(fname);
    num_images = numel(info);
    
    % to obtain the timing of the frames
    frames = voltageRec(:,3);
    meanFrame = nanmean(frames);
    frames(frames<meanFrame) = 0;
    frames(frames>meanFrame) = 1;
    [~,locsFrames] = findpeaks(frames);
    
    %to obtain the activity of the neurons
    unitVals = zeros(numNeuron, num_images);
    for k = 1:num_images
        Im = imread(fname, k);
        for nn = 1:numNeuron
            auxMask = mask;
            auxMask(auxMask~=nn) = 0;
            posx = find(sum(auxMask,1)~=0);
            posy = find(sum(auxMask,2)~=0);
            neuronMask = auxMask(posy(1):posy(end), posx(1):posx(end));
            Imd = double(Im(posy(1):posy(end),posx(1):posx(end)));
            unitVals(nn, k) = nansum(nansum(Imd.* neuronMask/nansum(neuronMask)));
        end
    end
    for nn = 1:numNeuron
        auxIndex = find(locs(nn)<=locsFrames, 1, 'first');
        minIndex = squeeze(max([auxIndex - wd, 1]));
        maxIndex = squeeze(min([auxIndex + wd, size(locsFrames,1) - 1, size(unitVals,2)]));
        subplot(numSubp,numSubp,nn)
        plot(locsFrames(minIndex:maxIndex), unitVals(nn, minIndex:maxIndex))
        vline(locsFrames(auxIndex), 'r-')
        title(['neuron: ', int2str(nn)])
    end
end