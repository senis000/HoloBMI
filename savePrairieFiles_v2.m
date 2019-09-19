function savePrairieFiles_v2(path_data, pl, expt_str)
% function to set the paths for saving the holo_stim_test
    %path_data.im - base directory where imaging data goes.
    savePathPrairieBase = fullfile(path_data.im, expt_str); 
    if ~exist(savePathPrairieBase, 'dir')
        mkdir(savePathPrairieBase);
    end
    saveCommand = "-p " + savePathPrairieBase;
    pl.SendScriptCommands(saveCommand);
    saveFilePrairie = "-fn Tseries " + expt_str + "_" + datestr(datetime('now'), 'yymmddTHHMMSS');
    pl.SendScriptCommands(saveFilePrairie);
    % run tseries with holoStim one neuron at a time ~2min and save the images.
end