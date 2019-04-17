function BMIAcqnvsPrairie(animal, day, E1, E2, T1, frameRate)
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
    ensemble = [E1, E2];
    numberNeurons = length(ensemble); 
    ensThresh = 0.95; % essentially top fraction (%) of what each cell can contribute to the BMI
    relaxationTime = 5;  % there can't be another hit in this many sec
    durationTrial = 30; % maximum time (in sec) that mice have for a trial
    movingAverage = 1; % Moving average of frames to calculate BMI (in sec)
    timeout = 5; %seconds of timeout if no hit in duration trial (sec)
    expectedLengthExperiment = 1*60*60*frameRate; % in frames
    baseLength = 2*60; % Period at the begginig without BMI to establish BL 

    savePath = ["F:/VivekNuria/", animal, "/",  day, "/"];
    
    % values of parameters in frames
    baseFrames = round(baseLength * frameRate);
    movingAverageFrames = round(movingAverage * frameRate);
    relaxationFrames = round(relaxationTime * frameRate);
    timeoutFrames = round(timeout * frameRate);

    %prairie view parameters
    chanIdx = 2; % green channel
    envPath = "F:/VivekNuria/"  ;%TODO set environment 

    % VTA parameters
    syncTime = 0.001; % duration of the TTL

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    global cursor miss hits trialEnd trialStart bmiAct
    

    
    %pre-allocating arrays
    expHistory = single(nan(numberNeurons, movingAverageFrames));  %define a windows buffer
    cursor = double(nan(1,expectedLengthExperiment));  %define a very long vector for cursor
    hits = single(zeros(1,expectedLengthExperiment));  %define a very long vector for hits
    miss = single(zeros(1,expectedLengthExperiment));  %define a very long vector for miss
    trialStart = single(zeros(1,expectedLengthExperiment));  %define a very long vector trialStart
    trialEnd = single(zeros(1,expectedLengthExperiment));  %define a very long vector trialEnd
    bmiAct = double(nan(numberNeurons, expectedLengthExperiment));


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

    tim = tic;

    % set the environment for the Time Series in PrairieView
    tslCommand = "tsl " + envPath;
    pl.SendScriptCommands(tslCommand);  %TODO check if this works

    
    %% Load Baseline variables
    % load onacid masks
    load([savePath, 'redcomp.mat']);
    neuronMask = reshape(full(AComp(:, ensemble)), px, py, numberNeurons);
    
    %% Create the file where to store info in case matlab crashes
    fileName = [savePath, 'bmiExp.dat'];
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    fwrite(fileID, bmiAct,'double');
    fwrite(fileID, cursor ,'double');
    fwrite(fileID, hits,'single');
    fwrite(fileID, miss,'single');
    fwrite(fileID, trialStart,'single');
    fwrite(fileID, trialEnd,'single');
    fwrite(fileID, T1,'single');
    fwrite(fileID, ensemble,'single');
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',{'double',size(bmiAct),'bmiAct'; ...
        'double',size(cursor),'cursor'; ...
        'single',size(hits),'hits'; ...
        'single',size(miss),'hits'; ...
        'single',size(trialStart),'trialStart'; ...
        'single',size(trialEnd),'trialEnd', 'repeat', 1); 
    m.Writable = true;
    
    %in case matlab crashes
    fileID = fopen([savePath, 'bmiExp.txt'],'wt');
    fprintf(fileID,'Length %d\n',round(expectedLengthExperiment));
    fprintf(fileID,'\n%6s %12s\r\n','E1','E2');
    A = [E1; E2];
    fprintf(fileID,'%6d %12d\r\n',A);
    fclose(fileID);
    
    %initialize the values of the memmap
    m.Data.trialStart = trialStart;
    m.Data.trialEnd = trialEnd;
    m.Data.hits = hits;
    m.Data.miss = miss;
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

            % Synchronization
            outputSingleScan(s,1);
            pause(syncTime)
            outputSingleScan(s,0);

            unitVals = obtainRoi(Im, neuronMask, com, ensemble); % function to obtain Rois values
            bmiAct(:,frame) = unitVals;
            m.Data.bmiAct(:,frame) = unitVals; % 1 ms  store info

            % update buffer of activity history
            expHistory(:, 1: end-1) = expHistory(:, 2:end);
            expHistory(:,end) = unitVals; 

            if counter == 0
                % Is it a new trial?
                if trialFlag && ~backtobaselineFlag
                    trialStart(frame) = 1;
                    m.Data.trialStart(frame) = 1;
                    trialHistory = trialHistory + 1;
                    trialFlag = 0;
                    %start running the timer again
                    tim = tic;
                    disp('New Trial!')
                end
                % calculate baseline activity and actual activity for the DFF

                signal = single(nanmean(expHistory, 2));
                if frame == 1
                    baseval = single(ones(numberNeurons,1)).*unitVals;
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
                m.Data.cursor(frame) = cursor(frame); % saving in memmap

                if backtobaselineFlag 
                    if cursor(frame) >= 1/2*T1   %TODO discuss!!!
                        backtobaselineFlag = 0;
                    end
                    tim = tic;  % to avoid false timeouts while it goes back to baseline
                else
                    if cursor(frame) <= T1      %if it hit the target
                        % VTA STIM
                        % Arduino pulse TODO
                        % a.writeDigitalPin("D6", 1); pause(0.010);a.writeDigitalPin("D6",0);

                        % update rewardHistory
                        rewardHistory = rewardHistory + 1;
                        disp(['Trial: ', num2str(trialHistory), 'Rewards: ', num2str(rewardHistory)]);
                        % update trials and hits vector
                        trialEnd(frame) = 1;
                        m.Data.trialEnd(frame) = 1;
                        hits(frame) = 1;
                        m.Data.hits(frame) = 1;

                        trialFlag = 1;
                        counter = relaxationFrames;
                        backtobaselineFlag = 1;

                    else
                        % update the tone to the new cursor
                        if toc(tim) > durationTrial
                            disp('Timeout')
                            trialEnd(frame) = 1;
                            m.Data.trialEnd(frame) = 1;
                            miss(frame) = 1;
                            m.Data.miss(frame) = 1;
                            trialFlag = 1;
                            counter = timeoutFrames;
                        end
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
    global cursor miss hits trialEnd trialStart bmiAct
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(savePath + "BMI_online.mat", 'cursor', 'miss', 'hits', 'trialEnd', ...
    'trialStart', 'bmiAct', 'animal', 'day', 'E1', 'E2', 'T1')
end




