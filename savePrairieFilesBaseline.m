function savePrairieFilesBaseline(savePath, pl)
% function to set the paths for saving the holo_stim_test

    savePathPrairie = fullfile(savePath, "im");
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    savePathPrairieBase = fullfile(savePathPrairie, "Baseline"); 
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    saveCommand = "-p " + savePathPrairieBase;
    pl.SendScriptCommands(saveCommand);
    saveFilePrairie = "-fn Tseries " + "Baseline_" + datestr(datetime('now'), 'yymmddTHHMMSS');
    pl.SendScriptCommands(saveFilePrairie);
    % run tseries with holoStim one neuron at a time ~2min and save the images.
end