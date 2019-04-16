function BaselineAcqnvsPrairie(animal, day, frameRate)
 %{
Function to acquire the baseline in a prairie scope
animal -> animal for the experiment
day -> day for the experiment
neuronMask -> matrix for spatial filters with px*py*unit 

%}

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    expectedLengthExperiment = 15*60*frameRate; % in frames
    
    savePath = ["F:/VivekNuria/", animal, "/",  day, "/"];

    %prairie view parameters
    chanIdx = 2; % green channel
    envPath = "F:/VivekNuria/"  ;%TODO set environment 

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    %pre-allocating arrays
    tempArray = single(nan(1,expectedLengthExperiment)); 

    %initializing flags and counters
    frame = 1; % initialize frames

    %% Prepare the nidaq
    s = daq.createSession('ni');
    addDigitalChannel(s,'dev5','Port0/Line0:0','OutputOnly');


    %% Prepare for Prairie
    % connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect()

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    % Prairie commands
    pl.SendScriptCommands("-srd True 0");
    pl.SendScriptCommands("-lbs True 0");

    lastFrame = zeros(px, py); % to compare with new incoming frames

    % set the environment for the Time Series in PrairieView
    tslCommand = "tsl " + envPath;
    pl.SendScriptCommands(tslCommand);  %TODO check if this works
    
    %% Load Baseline variables
    % load onacid masks
    load([savePath, 'redcomp.mat']);
    numberNeurons = size(AComp,2);
    neuronMask = reshape(full(AComp), px, py, numberNeurons);
    
    %% Create the file where to store the baseline
    baseActivity = zeros(numberNeurons, expectedLengthExperiment) + nan;
    fileName = [savePath, 'baselineActivity.dat'];
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    fwrite(fileID, baseActivity,'double');
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',{'double',size(baseActivity),'baseAct'}, 'repeat', 1); 
    m.Writable = true;
    %%
    %************************************************************************
    %*************************** RUN ********************************
    %************************************************************************

    %start the time_series scan
    pl.SendScriptCommands("-ts");   % TODO check if this actually works 


    while frame <= expectedLengthExperiment
        Im = pl.GetImage_2(chanIdx, px, py);
        if Im ~= lastFrame   
            lastFrame = Im;   % comparison and assignment takes ~4ms
            tempArray(frame) = datenum(clock); %TODO do we need it? do we want it?

            % Synchronization
            outputSingleScan(s,1);
            pause(syncTime)
            outputSingleScan(s,0);

            unitVals = obtainRoi(Im, neuronMask, com); % function to obtain Rois values 
            m.Data.baseAct(:,1) = unitVals; % 1 ms
            frame = frame + 1;
        end

    end
    pl.Disconnect();

end



