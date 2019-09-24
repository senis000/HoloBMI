function createBot_v2(save_file, posx,posy,r)
%{
Function to create a gpl file to be uploaded to the prairie view for BoT
save_file --> file where to save the .cfg files
posx,posy  -> position of the neurons
%}
    %% Parameters
    width = 8;  %half the width
    height = 8; %half the height
    
    %% print the first part of the text
    fileID = fopen(save_file,'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVBOTs>\n');
    
    % print for each point
    formatSpec = ['  <Region x="%d" y="%d" width="%d" height="%d" ', ...
        'angle="0" channel="2" shape="1" measurement="0" />\n'];
    
    for ind = 1:length(posx)
        %TODO change from pixelspace to motorspace
%         content = [posx(ind)-width, posy(ind)-height, width*2, height*2];
        content = [posx(ind)-r(ind), posy(ind)-r(ind), r(ind)*2, r(ind)*2];
        fprintf(fileID, formatSpec,content);
    end
    fprintf(fileID,'</PVBOTs>\n');
    fclose(fileID);
end