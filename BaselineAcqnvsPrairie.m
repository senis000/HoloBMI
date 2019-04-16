function BaselineAcqnvsPrairie(animal, day, neuronMask, frameRate)
 %{
Function to acquire the baseline in a prairie scope
animal -> animal for the experiment
day -> day for the experiment
neuronMask -> matrix for spatial filters with px*py*unit 

TODO we want to remove all neurons but the ensemble neurons from the
neuronMask

%}

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    expectedLengthExperiment = 15*60*frameRate; % in frames
    
    savePath = "F:/VivekNuria" + animal + day;

    %prairie view parameters
    chanIdx = 2; % green channel
    envPath = "F:/VivekNuria/"  ;%TODO set environment 

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************

    %pre-allocating arrays
    expHistory = single(nan(units, movingAverageFrames));  %define a windows buffer
    cursor = single(nan(1,expectedLengthExperiment));  %define a very long vector for cursor
    tempArray = single(nan(1,expectedLengthExperiment)); 
    miss = [];
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

            unitVals = obtainRoi(units, Im, neuronMask); % function to obtain Rois values

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
                    tim = tic;
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
                    if cursor(frame) >= 1/2*T1
                        backtobaselineFlag = 0;
                    end
                    tim = tic;  % to avoid false timeouts while it goes back to baseline
                else
                    if cursor(frame) <= T1 && ~motionFlag      %if it hit the target
                        % VTA STIM
                        %PYCONTROL PULSE!!!
    %                     outputSingleScan(s,1);
    %                     pause(VTATime)  %TODO better way to do this?
    %                     outputSingleScan(s,0);
                        % update rewardHistory
                        rewardHistory = rewardHistory + 1;
                        disp(['Trial: ', num2str(trialHistory), 'Rewards: ', num2str(rewardHistory)]);
                        % update trials and hits vector
                        trialEnd(end+1) = frame;
                        hits(end+1) = frame;
                        trialFlag = 1;
                        counter = relaxationFrames;
                        backtobaselineFlag = 1;

                    else
                        % update the tone to the new cursor
                        if toc(tim) > durationTrial
                            disp('Timeout')
                            trialEnd(end+1) = frame;
                            miss(end+1) = frame;
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
    global cursor miss hits trialEnd trialStart expHistory
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(savePath + "BMI_online.mat", 'cursor', 'miss', 'hits', 'trialEnd', ...
    'trialStart', 'expHistory', 'animal', 'day', 'E1', 'E2', 'T1')
end

function unitVals=obtainRoi(units, Im, neuronMask)
    for u = 1:units
        unitVals(u) = nanmean(nanmean(Im.*neuronMask(u,:,:)));
    end
end


