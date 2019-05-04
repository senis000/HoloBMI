function plotHoloStimTimeLock(botData, voltageRec, wd)

if nargin < 3
    wd = 100;
end

holoStim = voltageRec(:,8);
meanHoloStim = nanmean(holoStim);
holoStim(holoStim<meanHoloStim) = 0;
holoStim(holoStim>meanHoloStim) = 1;
[~,locs] = findpeaks(holoStim);

% Todo, check if this is needed or just a temporal delay in the recoding
temp = (botData(:,1) - botData(1,1))*1000;  %botData(:,1);
figure()
numSubp = ceil(sqrt(size(botData,2)-1));
for nn = 2:(size(botData,2))
    auxIndex = find(locs(nn)<=temp, 1, 'first');
    minIndex = max(auxIndex - wd, 1);
    maxIndex = min(auxIndex + wd, size(botData,1));
    subplot(numSubp,numSubp,nn-1)
    plot(temp(minIndex:maxIndex), botData(minIndex:maxIndex, nn))
    vline(temp(auxIndex), 'r-')
    title(['neuron: ', int2str(nn)])
       
end