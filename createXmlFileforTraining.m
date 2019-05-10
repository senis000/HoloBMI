function createXmlFileforTraining(savePath, power, IHSImean, IHSIrange, totalNumStims, expectedLengthExperiment, varName)
%{
Function to create a xml file to be uploaded to the prairie view
savePAth --> path where to save the .xml file
interholostim interval (IHSI) -> average time between stims
totalNumStims -> total number of stims
range -> dynamic range of IHSI
power 0.4 is a 100 in the gui
%}


    %% Parameters of the GUI
    Iter = 1;
    IterDelay = 1000.00;
    UncagingLaser = "Monaco";
    UncagingLaserPower = power;
    InterPointDelay =  0.12;
    Duration = 5;
    SpiralRevolutions = 15; 
    Repetitions = 8;  
    AllPointsAtOnce = "False";
    group = 1;
    numberNeuronsGroup = 4;
    
    
    %% Prepare vector of stims delay
    [~, delays] = createVectorHolo(1, expectedLengthExperiment, IHSImean, IHSIrange, false);
    delays = round(delays*1000);
    if length(delays) > totalNumStims
        delays = delays(1:totalNumStims);
    else
        disp('Not enough time to do all stims!!! ')
    end
    
    %% prepare file

    % print the first part of the text
    fileID = fopen(fullfile(savePath, [varName, 'groupholostim.xml']),'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVSavedMarkPointSeriesElements Iterations="%d" IterationDelay="%.2f">\n', Iter, IterDelay);
    
    % print for each point
    formatSpec1 = ['  <PVMarkPointElement Repetitions="%d" UncagingLaser="%s" UncagingLaserPower="%.1f"', ...
        ' TriggerFrequency="None" TriggerSelection="None" TriggerCount="1"', ...
        ' AsyncSyncFrequency="None" VoltageOutputCategoryName="None" VoltageRecCategoryName="None"', ...
        ' parameterSet="CurrentSettings">\n'];
    formatSpec2 = ['    <PVGalvoPointElement InitialDelay="%.2f"', ...
        ' InterPointDelay="%d" Duration="%d" SpiralRevolutions="%d" AllPointsAtOnce="%s"', ...
        ' Points="Group %d" Indices="1-%d" />\n'];
    
    for ind = 1:totalNumStims
        %TODO change from pixelspace to motorspace
        content1 = [delays(ind), InterPointDelay, Duration, SpiralRevolutions]; 
        content2 = [group, numberNeuronsGroup];
        fprintf(fileID, formatSpec1, Repetitions, UncagingLaser, UncagingLaserPower);
        fprintf(fileID, formatSpec2, content1, AllPointsAtOnce, content2);
        fprintf(fileID, '  </PVMarkPointElement>\n');
    end
    fprintf(fileID,'</PVSavedMarkPointSeriesElements>\n');
    fclose(fileID);
end