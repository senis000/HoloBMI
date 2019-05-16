function BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMask, onacid_bool, frameRate)
 %{
Function to acquire the baseline in a prairie scope
animal -> animal for the experiment
day -> day for the experiment
neuronMask -> matrix for spatial filters with px*py*unit 

%}

    if nargin<6
        onacid_bool = false;
    end

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    expectedLengthExperiment = 27000 * 3%ceil(3*60*frameRate); % in frames
    
    savePath = fullfile(folder, animal,  day);
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    %prairie view parameters
    chanIdx = 2; % green channel

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    global pl baseActivity
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(savePath));  %in case of ctrl-c it will lunch cleanmeup

    %% Prepare the nidaq
    s = daq.createSession('ni');
    addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
    ni_out = [0 0 0]; 
    outputSingleScan(s,ni_out);%set   
    ni_getimage = [1 0 0]; 


    %% Prepare for Prairie
    % connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect()
    
    % pause needed for prairie to respond
    pause(1)

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    % Prairie commands
    pl.SendScriptCommands("-srd True 0");
    pl.SendScriptCommands("-lbs True 0");
    
    %define where to save the file
    savePathPrairie = fullfile(savePath, "im");
    if ~exist(savePathPrairie, 'dir')
        mkdir(savePathPrairie);
    end
    savePrairieFiles(savePath, pl, "baseline")

    lastFrame = zeros(px, py); % to compare with new incoming frames

    % set the environment for the Time Series in PrairieView
    loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_15.env");
    pl.SendScriptCommands(loadCommand);  
    
    %% Load Baseline variables

    % create smaller versions of the spatial filter
    if onacid_bool
        numberNeurons = size(AComp,2);
        strcMask = obtainStrcMask(AComp, px, py);
    else
        numberNeurons = max(max(holoMask));
        strcMask = obtainStrcMaskfromMask(holoMask);
    end
    
    
    %% Create the file where to store the baseline
    baseActivity = zeros(numberNeurons, expectedLengthExperiment) + nan;
    fileName = fullfile(savePath, 'baselineActivity.dat');
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    if ~exist(fileName, 'file')
        disp('file does not exist. Memmap will not be saved')
    end
    fwrite(fileID, baseActivity,'double');
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',{'double',size(baseActivity),'baseAct'}, 'repeat', 1); 
    m.Writable = true;
    %%
    %************************************************************************
    %*************************** RUN ********************************
    %************************************************************************
    frame = 1; % initialize frames
    %start the time_series scan
    pl.SendScriptCommands("-ts");  
    disp('sent -ts'); 
    pause(2);  %empirically discovered time for the prairie to start gears
    counterSame = 0;
    disp('Starting baseline acquisition')
    while counterSame < 500
        Im = pl.GetImage_2(chanIdx, px, py);
        if ~isequal(Im,lastFrame)   
            tic
            lastFrame = Im;   % comparison and assignment takes ~4ms
            outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,[0 0 0]);
           
            unitVals = obtainRoi(Im, strcMask); % function to obtain Rois values 
            baseActivity(:,frame) = unitVals;
            m.Data.baseAct(:,frame) = unitVals; % 1 ms
            frame = frame + 1;
            counterSame = 0;
            t = toc;
            if t < 1/(frameRate*1.2)
                pause(1/(frameRate*1.2) -t)
            end
        else
            counterSame = counterSame + 1;
        end

    end
   % pl.Disconnect();

end

function cleanMeUp(savePath)
    global pl baseActivity
    disp('cleaning')
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(fullfile(savePath, "BaselineOnline" + datestr(datetime('now'), 'yymmddTHHMMSS') + ".mat"), 'baseActivity')
    if pl.Connected()
        pl.Disconnect();
    end
end



