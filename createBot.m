function createBot(savePath, posx,posy)
%{
Function to create a gpl file to be uploaded to the prairie view for BoT
savePath --> path where to save the .gpl file
posx,posy  -> position of the neurons
%}
    %% Parameters
    width = 10;
    height = 10;

    %% print the first part of the text
    fileID = fopen(fullfile(savePath, 'Bot.gpl'),'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVBOTs>\n');
    
    % print for each point
    formatSpec = ['  <Region x="%d" y="%d" width="%d" height="%d" ', ...
        'angle="0" channel="2" shape="1" measurement="0" />\n'];
    
    for ind = 1:length(posx)
        %TODO change from pixelspace to motorspace
        content = [posx(ind), posy(ind), width, height];
        fprintf(fileID, formatSpec,content);
    end
    fprintf(fileID,'</PVBOTs>\n');
    fclose(fileID);
end