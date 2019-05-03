function createGplFile(savePath, holoMask, posz, pixelSize, varName)
%{
Function to create a gpl file to be uploaded to the prairie view
savePAth --> path where to save the .gpl file
holoMask -> mask of the red neurons to be activated
%}

    if nargin <5
        varName = '';
    end

    %% Parameters
    UncagingLaserPower = 0;
    Duration = 100;
    SpiralSize = 0.2;
    SpiralRevolutions = 5; 
    
    conversionValue = 5.06666666666667;
    
    %% prepare file
    % obtain positions in the mask
    [posx, posy] = findCenter(holoMask);
    
    % change from pixel space to motor space
    posx = 2*conversionValue/pixelSize*posx - conversionValue;
    posy = -2*conversionValue/pixelSize*posy + conversionValue;
    
    % print the first part of the text
    fileID = fopen(fullfile(savePath, [varName, 'holoMask.gpl']),'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVGalvoPointList>\n');
    
    % print for each point
    formatSpec = ['  <PVGalvoPoint X="%f" Y="%f" Name="Point %d" Index="%d"', ...
        ' ActivityType="MarkPoints" UncagingLaser="None" UncagingLaserPower="%d"', ...
        ' Duration="%d" IsSpiral="True" SpiralSize="%f" SpiralRevolutions="%d" Z="%.2f" />\n'];
    
    for ind = 1:length(posx)
        %TODO change from pixelspace to motorspace
        content = [posx(ind), posy(ind), ind, ind-1, UncagingLaserPower, Duration, ...
            SpiralSize, SpiralRevolutions, posz];
        fprintf(fileID, formatSpec,content);
    end
    fprintf(fileID,'</PVGalvoPointList>\n');
    fclose(fileID);
end