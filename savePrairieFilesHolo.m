function savePrairieFilesHolo(savePath)
% function to set the paths for saving the holo_stim_test
%savePathPrairie = savePath + im
%savePathPrairieHolo = savePath + im + holostim
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    pause(2);
    savePathPrairie = fullfile(savePath, "im");
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    savePathPrairieHolo = fullfile(savePathPrairie, "holostim"); 
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    saveCommand = "-p " + savePathPrairieHolo;
    pl.SendScriptCommands(saveCommand);
    saveFilePrairie = "-fn Tseries " + "testholo_" + datestr(datetime('now'), 'yymmddTHHMMSS');
    pl.SendScriptCommands(saveFilePrairie);
    % run tseries with holoStim one neuron at a time ~2min and save the images.
    pl.Disconnect()
end