function createXplFile(savePath)
%{
Function to create a gpl file to be uploaded to the prairie view
savePAth --> path where to save the .gpl file
holoMask -> mask of the red neurons to be activated
%}

    %% Parameters
    Iterations = 1;
    IterationDelay = 1000.00;
    UncagingLaser = 'fidelityHP';
    UncagingLaserPower = 0.4;
    InitialDelay = 0.22;
    SpiralSize = 2.1067852798579;
    SpiralRevolutions = 5; 
    
    %% prepare file
    % obtain positions in the mask
    [posx, posy] = findCenter(holoMask);
    % print the first part of the text
    fileID = fopen([savePath, 'holoMask.gpl'],'wt');
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