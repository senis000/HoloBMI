function HoloAcqnvsPrairie_v2(path_data, expt_str, mask, expectedLengthExperiment)
%path_data.savePath
%path_data.im
% HoloAcqnvsPrairie_v2(folder, animal, day, mask)
 %{
Function to acquire the holostim in a prairie scope
animal -> animal for the experiment
day -> day for the experiment
neuronMask -> matrix for spatial filters with px*py*unit 

%}

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
%     expectedLengthExperiment = 3*7000%ceil(3*60*frameRate); % in frames
    %prairie view parameters
    chanIdx = 2; % green channel
%     envPath = fullfile(folder,"utils", "Tseries_VivekNuria_15.env")  ;%TODO set environment 

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    global pl holoActivity
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(path_data.savePath, expt_str));  %in case of ctrl-c it will lunch cleanmeup

    %% Prepare the nidaq
%     s = daq.createSession('ni');
%     addDigitalChannel(s,'dev5','Port0/Line0:0','OutputOnly');


    %% Prepare for Prairie
    % connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect()
    
    % pause needed for prairie to respond
    pause(2)

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    % Prairie commands
    pl.SendScriptCommands("-srd True 0");
    pl.SendScriptCommands("-lbs True 0");
    
    %define where to save the file
%     savePathPrairie = path_data.im;
    savePrairieFiles_v2(path_data, pl, expt_str)

    lastFrame = zeros(px, py); % to compare with new incoming frames

    % set the environment for the Time Series in PrairieView
%     tslCommand = "-tsl " + envPath;
%     pl.SendScriptCommands(tslCommand);  
    
    %% Load variables

    numberNeurons = max(max(mask));
    % create smaller versions of the spatial filter
    strcMask = obtainStrcMaskfromMask(mask);
    
    %% Prepare the nidaq
    clear s
    s = daq.createSession('ni');
    addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
    ni_out = [0 0 0]; 
    outputSingleScan(s,ni_out);%set   
    ni_getimage = [1 0 0]; 
    
    %% Create the file where to store the holostim
    holoActivity = zeros(numberNeurons, expectedLengthExperiment) + nan;
    fileName = fullfile(savePath, [expt_str '.dat']);
%     fileName = fullfile(savePath, 'holoActivity.dat');
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    if ~exist(fileName, 'file')
        disp('file does not exist. Memmap will not be saved')
    end
    fwrite(fileID, holoActivity,'double');
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',{'double',size(holoActivity),'holoAct'}, 'repeat', 1); 
    m.Writable = true;
    %%
    %************************************************************************
    %*************************** RUN ********************************
    %************************************************************************
    frame = 1; % initialize frames
    %start the time_series scan
    pl.SendScriptCommands("-ts"); 
    pause(1);  %empirically discovered time for the prairie to start gears

    disp('Starting holo acquisition')
    counterSame = 0;
    while counterSame < 500
        Im = pl.GetImage_2(chanIdx, px, py);
        if ~isequal(Im,lastFrame) 
            tic;
            lastFrame = Im;   % comparison and assignment takes ~4ms
            outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,[0 0 0]);
            
            unitVals = obtainRoi(Im, strcMask); % function to obtain Rois values 
            holoActivity(:,frame) = unitVals;
            m.Data.holoAct(:,frame) = unitVals; % 1 ms
            frame = frame + 1;
            counterSame = 0;
            t = toc;
            if t< 0.02
                pause(0.02-t)
            end
        else
            counterSame = counterSame + 1;
        end

    end
   % pl.Disconnect();

end

function cleanMeUp(savePath, expt_str)
    global pl holoActivity
    disp('cleaning')
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
%     save(fullfile(savePath, "holoOnline" + datestr(datetime('now'), 'yymmddTHHMMSS') + ".mat"), 'holoActivity')
    save(fullfile(savePath, expt_str + datestr(datetime('now'), 'yymmddTHHMMSS') + ".mat"), 'holoActivity')
    if pl.Connected()
        pl.Disconnect();
    end
end



