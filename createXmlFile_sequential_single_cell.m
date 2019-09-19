function createXmlFile_sequential_single_cell(xml_path, params, stim_sequence)

%{
Function to create a xml file to be uploaded to the prairie view
xml_path --> file's full path to be saved to.
savePAth --> path where to save the .xml file
numberNeurons -> How many neurons to randomly activate
reps --> amount of stim per neuron
flagRandom --> adding randomness to the targets of stimulation
Random only happens within repetitions
power 0.4 is a 100 in the gui
%}
    p = params; %rename for brevity
    
    %% prepare file

    % print the first part of the text
    fileID = fopen(xml_path,'wt');
    fprintf(fileID,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fileID,'<PVSavedMarkPointSeriesElements Iterations="%d" IterationDelay="%.2f">\n', p.Iter, p.IterDelay);
    
    % print for each point
    formatSpec1 = ['  <PVMarkPointElement Repetitions="%d" UncagingLaser="%s" UncagingLaserPower="%.3f"', ...
        ' TriggerFrequency="None" TriggerSelection="None" TriggerCount="1"', ...
        ' AsyncSyncFrequency="None" VoltageOutputCategoryName="None" VoltageRecCategoryName="None"', ...
        ' parameterSet="CurrentSettings">\n'];
    formatSpec2 = ['    <PVGalvoPointElement InitialDelay="%.2f"', ...
        ' InterPointDelay="%d" Duration="%d" SpiralRevolutions="%d" AllPointsAtOnce="%s"', ...
        ' Points="Point %d" Indices="%d" />\n'];
    
    for i_stim = 1:length(stim_sequence)
        n = stim_sequence(i_stim); %which neuron to stim
        %TODO change from pixelspace to motorspace
        fprintf(fileID, formatSpec1, p.RepetitionsVector(n), p.UncagingLaser, p.PowerVector(n));
        content1 = [p.InitialDelayVector(n), p.InterPointDelayVector(n), p.DurationVector(n), p.SpiralVector(n)]; 
        content2 = [n, n];        
        fprintf(fileID, formatSpec2, content1, p.AllPointsAtOnce, content2);
        fprintf(fileID, '  </PVMarkPointElement>\n');
    end
    fprintf(fileID,'</PVSavedMarkPointSeriesElements>\n');
    fclose(fileID);
end