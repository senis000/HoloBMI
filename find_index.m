function [uniqueNmat, uniqueNcaim, correlations] =find_index(com, CComp, YrA, red, holoActivity, voltageRec, auxtol)

    if nargin < 7
        auxtol = 4;
    end
    
    % initialize vars
    
	% find frames
    matframes = voltageRec(:,9);
    meanFrame = nanmean(matframes);
    matframes(matframes<meanFrame) = 0;
    matframes(matframes>meanFrame) = 1;
    [~,locsFramesMat] = findpeaks(matframes);
    
    % to obtain the timing of the frames
    frames = voltageRec(:,3);
    meanFrame = nanmean(frames);
    frames(frames<meanFrame) = 0;
    frames(frames>meanFrame) = 1;
    [~,locsFramesIm] = findpeaks(frames);
    if length(locsFramesIm) > size(CComp,2)
        locsFramesIm = locsFramesIm(1:size(CComp,2));
    end
    
    % find distances
    dist = pdist2(red', com);
    
    dist(dist>auxtol) = nan;
    [nmat, ncaim] = find(~isnan(dist));
    
    uniqueNmat = unique(nmat);
    uniqueNcaim = zeros(1, size(uniqueNmat,1));
    correlations = zeros(1, size(uniqueNmat,1));
    CComp = CComp + YrA;
    for indm = 1:length(uniqueNmat)
        possible_match = ncaim(find(nmat==uniqueNmat(indm)));
        neurcor = zeros(length(possible_match), 1);
        for indc = 1: length(possible_match)
            resampled = sort(unique([locsFramesIm; locsFramesMat]));
            auxH = holoActivity(uniqueNmat(indm),1:size(locsFramesMat,1));
            auxH = (auxH - nanmean(auxH))./nanmean(auxH);
            auxC = CComp(possible_match(indc),:)/10000;
            vC = interp1(locsFramesIm,auxC,resampled);
            vH = interp1(locsFramesMat,auxH,resampled);
            cc = corrcoef(vH, vC, 'Rows','complete');
            neurcor(indc) = cc(1,2);
        end
        [corval,maxcor] = nanmax(neurcor);
        uniqueNcaim(indm) = possible_match(maxcor);
        correlations(indm) = corval;
        fprintf('Neuron %d of matlab corresponds to neuron %d of caiman with a correlation %.2f \n', uniqueNmat(indm), uniqueNcaim(indm), corval);
    end
            
