function [save_file] = createGplFile_v2(savePath, params, posx, posy, posz_scalar, r, num_pixels, zoom, varName)
%{
Function to create a gpl file to be uploaded to the prairie view
savePath --> path where to save the .gpl file
posx --> vector of x coordinates
poxy --> vector of y coordinates
r --> vector of roi size
%}

    if nargin <9
        varName = '';
    end

    %% Parameters
%     UncagingLaserPower = 0.4;
%     Duration = 100;
%     SpiralSize = 0.3;
%     SpiralRevolutions = 5; 
    if(zoom == 1.5)
        conversionValue = 5.06666666666667;
    elseif(zoom == 2)
        conversionValue = 3.8;
        disp('zoom2');
        %512 > 3.8
        %0 -> -3.8
    else
        print('please calibrate the zoom')
    end
    
    %% prepare file
    
    % change from pixel space to motor space
    posx = 2*conversionValue/num_pixels*posx - conversionValue;
    posy = -2*conversionValue/num_pixels*posy + conversionValue;
    
    % print the first part of the text
    save_file = fullfile(savePath, [varName, 'holoMask.gpl']);
    fileID = fopen(save_file,'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVGalvoPointList>\n');
    
    % print for each point
    formatSpec = ['  <PVGalvoPoint X="%f" Y="%f" Name="Point %d" Index="%d"', ...
        ' ActivityType="MarkPoints" UncagingLaser="None" UncagingLaserPower="%d"', ...
        ' Duration="%d" IsSpiral="True" SpiralSize="%f" SpiralRevolutions="%d" Z="%.2f" />\n'];
    
    num_roi = length(posx); 
    for ind = 1:num_roi
        %TODO change from pixelspace to motorspace
        content = [posx(ind), posy(ind), ind, ind-1, params(ind).UncagingLaserPower, params(ind).Duration, ...
            params(ind).SpiralSize, params(ind).SpiralRevolutions, posz_scalar];
        fprintf(fileID, formatSpec,content);
    end
    fprintf(fileID,'</PVGalvoPointList>\n');
    fclose(fileID);
end