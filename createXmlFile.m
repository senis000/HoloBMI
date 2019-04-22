function createXmlFile(savePath, numberNeurons, reps, flagRandom)
%{
Function to create a gpl file to be uploaded to the prairie view
savePAth --> path where to save the .gpl file
numberNeurons -> How many neurons to randomly activate
reps --> amount of stim per neuron
flagRandom --> adding randomness to the targets of stimulation
Random only happens within repetitions
%}

    if nargin < 4
        flagRandom = True;
    end

    %% Parameters of the GUI
    Iter = 1;
    IterDelay = 1000.00;
    UncagingLaser = "fidelityHP";
    UncagingLaserPower = 0.4;
    InitialDelay = 0.22;
    InterPointDelay = 4800;
    Duration = 200;
    SpiralRevolutions = 5; 
    Repetitions = 1;  
    AllPointsAtOnce = "False";
    
    %% Array of random holostims 
    holovector = nan(numberNeurons,reps);
    for rr = 1:reps
        if flagRandom
            holovector(:,rr) = randsample(numberNeurons,numberNeurons);
        else
            holovector(:,rr) = 1:numberNeurons;
        end
    end
    holovector = holovector(:);
    
    
    %% prepare file

    % print the first part of the text
    fileID = fopen([savePath, 'holoMask.xml'],'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVSavedMarkPointSeriesElements Iterations="%d" IterationDelay="%.2f">\n', Iter, IterDelay);
    
    % print for each point
    formatSpec1 = ['  <PVMarkPointElement Repetitions="%d" UncagingLaser="%s" UncagingLaserPower="%.1f"', ...
        ' TriggerFrequency="None" TriggerSelection="None" TriggerCount="1"', ...
        ' AsyncSyncFrequency="None" VoltageOutputCategoryName="None" VoltageRecCategoryName="None"', ...
        ' parameterSet="CurrentSettings">\n'];
    formatSpec2 = ['    <PVGalvoPointElement InitialDelay="%.2f"', ...
        ' InterPointDelay="%d" Duration="%d" SpiralRevolutions="%d" AllPointsAtOnce="%s"', ...
        ' Points="Point %d" Indices="%d" />\n'];
    
    for ind = 1:numberNeurons*reps
        %TODO change from pixelspace to motorspace
        content1 = [InitialDelay, InterPointDelay, Duration, SpiralRevolutions]; 
        content2 = [holovector(ind), ind];
        fprintf(fileID, formatSpec1, Repetitions, UncagingLaser, UncagingLaserPower);
        fprintf(fileID, formatSpec2, content1, AllPointsAtOnce, content2);
        fprintf(fileID, '  </PVMarkPointElement>\n');
    end
    fprintf(fileID,'</PVSavedMarkPointSeriesElements>\n');
    fclose(fileID);
end