function savePrairieFiles(savePath, pl, expt_str)
% function to set the paths for saving the holo_stim_test

    savePathPrairie = fullfile(savePath, "im");
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    savePathPrairieBase = fullfile(savePathPrairie, expt_str); 
    if ~exist(savePathPrairieBase, 'dir')
        mkdir(savePathPrairieBase);
    end
    saveCommand = "-p " + savePathPrairieBase;
    pl.SendScriptCommands(saveCommand);
    saveFilePrairie = "-fn Tseries " + expt_str + "_" + datestr(datetime('now'), 'yymmddTHHMMSS');
    pl.SendScriptCommands(saveFilePrairie);
    % run tseries with holoStim one neuron at a time ~2min and save the images.
end