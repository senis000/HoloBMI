function BMIAcqnvsPrairienoTrials(animal, day, E1, E2, T1, frameRate)
    %{
    Function to acquire the BMI in a prairie scope
    animal -> animal for the experiment
    day -> day for the experiment
    neuronMask -> matrix for spatial filters with px*py*unit 
    and nan otherwise
    E2 = [1 2 3 4]; index in neuronMask for the ensembles. E2 being the one
    that has to increase
    E1 = [5 6 7 8]; 
    T1 = target for reward

%}


    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    % BMI parameters 
    %frameRate = 30; % TODO check if it can be obtained from prairie
    units = length(E1)+length(E2); 
    ensThresh = 0.75; % essentially top fraction (%) of what each cell can contribute to the BMI
    relaxationTime = 5;  % there can't be another hit in this many sec
    movingAverage = 1; % Moving average of frames to calculate BMI (in sec)
    expectedLengthExperiment = 1*60*60*frameRate; % in frames
    baseLength = 2*60; % Period at the begginig without BMI to establish BL 

    savePath = ["F:/VivekNuria/", animal, "/",  day, "/"];

    % values of parameters in frames
    baseFrames = round(baseLength * frameRate);
    movingAverageFrames = round(movingAverage * frameRate);
    relaxationFrames = round(relaxationTime * frameRate);

    %prairie view parameters
    chanIdx = 2; % green channel
    envPath = "F:/VivekNuria/"  ;%TODO set environment 

    % VTA parameters
    syncTime = 0.001; % duration of the TTL

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    global cursor hits trialEnd trialStart
    
    %pre-allocating arrays
    expHistory = single(nan(units, movingAverageFrames));  %define a windows buffer
    cursor = single(nan(1,expectedLengthExperiment));  %define a very long vector for cursor
    tempArray = single(nan(1,expectedLengthExperiment)); 
    hits = [];
    trialEnd = [];
    trialStart = [];


    %initializing flags and counters
    rewardHistory = 0; 
    trialHistory = 0;
    trialFlag = 1;
    counter = 40;
    backtobaselineFlag = 0;
    frame = 1; % initialize frames

    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(savePath, animal, day, E1, E2, T1));  %in case of ctrl-c it will lunch cleanmeup


    %% Prepare the nidaq
    s = daq.createSession('ni');
    addDigitalChannel(s,'dev5','Port0/Line0:0','OutputOnly');
    
    %% Prepare the arduino
    a = arduino('COM15', 'Uno');  
    %a.writeDigitalPin("D6", 1); pause(0.005);a.writeDigitalPin("D6",0); 
    %TODO check port, pin and pause time


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
            tempArray(frame) = datenum(clock);


            % Synchronization
            outputSingleScan(s,1);
            pause(syncTime)
            outputSingleScan(s,0);

            unitVals = obtainRoi(Im, neuronMask, com, [E1, E2]); % function to obtain Rois values

            % update buffer of activity history
            expHistory(:, 1: end-1) = expHistory(:, 2:end);
            expHistory(:,end) = unitVals; 

            if counter == 0
                % Is it a new trial?
                if trialFlag && ~backtobaselineFlag
                    trialStart(end+1) = frame;
                    trialHistory = trialHistory + 1;
                    trialFlag = 0;
                    %start running the timer again
                    disp('New Trial!')
                end
                % calculate baseline activity and actual activity for the DFF

                signal = single(nanmean(expHistory, 2));
                if frame == 1
                    baseval = single(ones(units,1)).*unitVals;
                elseif frame <= baseFrames
                    baseval = (baseval*(frame - 1) + signal)./frame;
                else
                    baseval = (baseval*(baseFrames - 1) + signal)./baseFrames;
                end

                % calculation of DFF
                dff = (signal - baseval) ./ baseval;
                dff(dff<T1*ensThresh) = T1*ensThresh; % limit the contribution of each neuron to a portion of the target
                % it is very unprobable (almost imposible) that one cell of E1 does
                % it on its own, buuut just in case:
                dff(dff>-T1*ensThresh) = -T1*ensThresh;

                cursor(frame) = nansum([nansum(dff(E1)),-nansum(dff(E2))]);

                if backtobaselineFlag 
                    if cursor(frame) >= 1/2*T1  %TODO discuss!!!
                        backtobaselineFlag = 0;
                    end
                else
                    if cursor(frame) <= T1      %if it hit the target
                        % VTA STIM
                        % Arduino pulse TODO
                        % a.writeDigitalPin("D6", 1); pause(0.010);a.writeDigitalPin("D6",0);
                        
                        % update rewardHistory
                        rewardHistory = rewardHistory + 1;
                        disp(['Trial: ', num2str(trialHistory), 'Rewards: ', num2str(rewardHistory)]);
                        % update trials and hits vector
                        trialEnd(end+1) = frame;
                        hits(end+1) = frame;
                        trialFlag = 1;
                        counter = relaxationFrames;
                        backtobaselineFlag = 1;
                    end
                end
            else        
                counter = counter - 1;
            end
            frame = frame + 1;
        end

    end
    pl.Disconnect();

end

% fires when main function terminates (normal, error or interruption)
function cleanMeUp(savePath, animal, day, E1, E2, T1)
    global cursor hits trialEnd trialStart 
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(savePath + "BMI_online.mat", 'cursor', 'hits', 'trialEnd', ...
    'trialStart', 'animal', 'day', 'E1', 'E2', 'T1')
end


