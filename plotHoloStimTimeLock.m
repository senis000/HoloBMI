function [ind, zVector] = plotHoloStimTimeLock(holoData, voltageRec, duration, wd)
% duration = 40
% wd = 1000

%{
Function to plot the response of neurons to the holo_stim
botData => data from the BOT file
voltageRec => data from the voltage recording
wd => window to plot before/after neuron
duration => including the repetitions
%}
%     if nargin < 4
%         wd = 100;
%     end
    
    preStim = 2*1000;
    postStim = 0.25*1000;

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
    
%     figure()
    numSubp = ceil(sqrt(size(holoData,1)-1));
    indSubp = min([size(holoData,1)-1, length(locs)]);
    meanVector = zeros(indSubp, 1);
    stdVector = zeros(indSubp, 1);
    zVector = zeros(indSubp, 1);
    
    holoData_dff = holoData;
    holoData_detrend = holoData_dff;
    for nn =1:indSubp
        holoData_dff(nn,:) = ...
            (holoData(nn,:)-nanmean(holoData(nn, :)))/nanmean(holoData(nn, :));
        
        %Detrend:
        d_i = holoData_dff(nn,:);
        %remove nans:
        d_i_valid = d_i(~isnan(d_i)); 
        d_i_detrend = detrend(d_i_valid);
        holoData_detrend(nn,:) = d_i;
        holoData_detrend(nn,~isnan(d_i)) = d_i_detrend;
    end
    
    z_data = {};
    for nn = 1: indSubp
        auxIndex = find(locs(nn)<=locsFrames, 1, 'first');
        minIndex = max(auxIndex - wd, 1);
        maxIndex = min(auxIndex + wd, size(locsFrames,1));
        minzIndex = max(auxIndex - preStim, 1);
        maxzIndex = min(auxIndex + postStim, size(holoData,2));
        meanVector(nn) = nanmean(holoData_detrend(nn, minzIndex:auxIndex),2);
        stdVector(nn) = nanstd(holoData_detrend(nn, :),0,2);
%         stdVector(nn) = nanstd(holoData_detrend(nn, minzIndex:auxIndex),0,2); %nanstd(holoData(nn, minzIndex:auxIndex),0,2);
        z_data{nn} = (holoData_detrend(nn, auxIndex:maxzIndex))/stdVector(nn); % - meanVector(nn
        zVector(nn) = nanmean(z_data{nn});
%         subplot(numSubp,numSubp,nn)
%         plot(locsFrames(minIndex:maxIndex)./1000, holoData_dff(nn, minIndex:maxIndex))
%         vline(locs(nn)/1000, 'r-')
%         title(['neuron: ', int2str(nn)])
    end
    
    [~, ind] = sort(zVector, 'descend');
    
    figure()
    for nn = 1: length(ind)
        auxIndex = find(locs(ind(nn))<=locsFrames, 1, 'first');
        minIndex = max(auxIndex - wd, 1);
        maxIndex = min(auxIndex + wd, size(locsFrames,1));
        subplot(numSubp,numSubp,nn)
        plot(locsFrames(minIndex:maxIndex)./1000, holoData_dff(ind(nn), minIndex:maxIndex))
        vline(locs(ind(nn))/1000, 'r-')
        title(['neuron: ', int2str(ind(nn))])
    end
    
% end

%%
% %65 vs 9
% nn=62;
% y_data = holoData_detrend(nn,:);
% auxIndex = find(locs(nn)<=locsFrames, 1, 'first')
% 
% h = figure; 
% plot(y_data); 
% vline(auxIndex)
% 
% h = figure;
% plot(z_data{nn})
% title(num2str(mean(z_data{nn}))); 
% vline(locs(nn)/1000, 'r-')

% %%
% h = figure;
% plot(holoData_detrend(nn,:)); 
% 
% %%
% h = figure;
% plot(holoData(nn,:)); 


% %%
% mean(z_data{nn})
% %%
% h = figure;
% plot(zVector(ind), '.-')
