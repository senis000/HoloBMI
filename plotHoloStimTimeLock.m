function [ind, zVector] = plotHoloStimTimeLock(holoData, voltageRec, duration, wd)
%{
Function to plot the response of neurons to the holo_stim
botData => data from the BOT file
voltageRec => data from the voltage recording
wd => window to plot before/after neuron
duration => including the repetitions
%}
    if nargin < 4
        wd = 100;
    end
    
    preStim = 5*1000;
    postStim = 2*1000;

    holoStim = voltageRec(:,8);
    meanHoloStim = nanmean(holoStim);
    holoStim(holoStim<meanHoloStim) = 0;
    holoStim(holoStim>meanHoloStim) = 1;
    [~,locs] = findpeaks(holoStim);
    
    if length(locs) > size(holoData,1)
        locs(find(diff(locs)<duration)+1) = [];
    end
    
    frames = voltageRec(:,9);
    meanFrame = nanmean(frames);
    frames(frames<meanFrame) = 0;
    frames(frames>meanFrame) = 1;
    [~,locsFrames] = findpeaks(frames);
    
    % Todo, check if this is needed or just a temporal delay in the recoding
    
    figure()
    numSubp = ceil(sqrt(size(holoData,1)-1));
    indSubp = min([size(holoData,1)-1, length(locs)]);
    meanVector = zeros(indSubp, 1);
    stdVector = zeros(indSubp, 1);
    zVector = zeros(indSubp, 1);
    
    for nn = 1: indSubp
        auxIndex = find(locs(nn)<=locsFrames, 1, 'first');
        minIndex = max(auxIndex - wd, 1);
        maxIndex = min(auxIndex + wd, size(holoData,2));
        minzIndex = max(auxIndex - preStim, 1);
        maxzIndex = min(auxIndex + postStim, size(holoData,2));
        meanVector(nn) = nanmean(holoData(nn, :),2);
        stdVector(nn) = nanstd(holoData(nn, :),0,2); %nanstd(holoData(nn, minzIndex:auxIndex),0,2);
        zVector(nn) = nanmean((holoData(nn, auxIndex:maxzIndex) - meanVector(nn))/stdVector(nn));
        subplot(numSubp,numSubp,nn)
        plot(round(locsFrames(minIndex:maxIndex)./100)./10, holoData(nn, minIndex:maxIndex)./nanmean(holoData(nn, :)))
        vline(locs(nn)/1000, 'r-')
        title(['neuron: ', int2str(nn)])
    end
    
    [~, ind] = sort(zVector, 'descend');
    
    figure()
    for nn = 1: length(ind)
        auxIndex = find(locs(ind(nn))<=locsFrames, 1, 'first');
        minIndex = max(auxIndex - wd, 1);
        maxIndex = min(auxIndex + wd, size(holoData,2));
        subplot(numSubp,numSubp,nn)
        plot(round(locsFrames(minIndex:maxIndex)./100)./10, holoData(ind(nn), minIndex:maxIndex)./nanmean(holoData(ind(nn), :)))
        vline(locs(ind(nn))/1000, 'r-')
        title(['neuron: ', int2str(ind(nn))])
    end
    
end
