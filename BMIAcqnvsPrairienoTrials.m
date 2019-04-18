function BMIAcqnvsPrairienoTrials(animal, day, E1Base, E2Base, T1, frameRate, vectorHolo)
    %{
    Function to acquire the BMI in a prairie scope
    animal -> animal for the experiment
    day -> day for the experiment
    neuronMask -> matrix for spatial filters with px*py*unit 
    and nan otherwise
    E2base = [10 2 13 24]; index in neuronMask for the ensembles. E2 being the one
    that has to increase
    E1Base = [54 6 17 28]; 
    They will change to a new 
    T1 = target for reward (positive)
    vectorHolo -> vector of scheduled holo stims
%}

    if nargin < 7
        vectorHolo = [];
    end

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    
    %% experiment FLAGS
    flagVTA = true;
    flagHolo = true;
    
    %% BMI parameters 
    %frameRate = 30; % TODO check if it can be obtained from prairie 
    relaxationTime = 0;  % there can't be another hit in this many sec
    back2Base = 1/2*T1; % cursor must be under this value to be able to hit again

    savePath = ['F:/VivekNuria/', animal, '/',  day, '/'];
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    % values of parameters in frames
    expectedLengthExperiment = 2*60*frameRate; % in frames
    baseFrames = round(2*60 * frameRate); % Period at the begginig without BMI to establish BL
    movingAverageFrames = 2;
    relaxationFrames = round(relaxationTime * frameRate);

    %% prairie view parameters
    chanIdx = 2; % green channel
    envPath = 'F:/VivekNuria/utils/'  ;%TODO set environment 

    %% VTA parameters
    shutterVTA = round(2*frameRate);
    syncVTA = 0.001; % duration of the TTL

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************
    
    global pl cursor hits trialStart bmiAct baseVector timeVector %TODO remove timeVector
    
    %obtain new E1, E2 and save AComp
    [E1, E2] = fromBaseline2Bmi(E1Base, E2Base);
    ensemble = [E1, E2];
    numberNeurons = length(ensemble);
    
    %pre-allocating arrays
    expHistory = single(nan(numberNeurons, movingAverageFrames));  %define a windows buffer
    cursor = double(nan(1,expectedLengthExperiment));  %define a very long vector for cursor
    hits = single(zeros(1,expectedLengthExperiment));  %define a very long vector for hits
    trialStart = single(zeros(1,expectedLengthExperiment));  %define a very long vector trialStart
    bmiAct = double(nan(numberNeurons, expectedLengthExperiment));
    baseVector = double(nan(numberNeurons,expectedLengthExperiment));  %define a very long vector for cursor

    %to debug!!! TODO REMOVE after debugging
    timeVector = double(nan(1,expectedLengthExperiment));  %define a very long vector for cursor

    %initializing general flags and counters 
    rewardHistory = 0; 
    trialHistory = 0;
    trialFlag = 1;
    nonBufferUpdateCounter = 40;  %counter when we dont want to update the buffer: 
    initFrameBase = nonBufferUpdateCounter + 1;
    %beginning of experiment and VTA stim
    BufferUpdateCounter = 0;
    
    backtobaselineFlag = 0;
    frame = 1; % initialize frames
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(savePath, animal, day, E1, E2, T1, E1Base, E2Base));  %in case of ctrl-c it will lunch cleanmeup

%     %% Prepare the nidaq
%     s = daq.createSession('ni');
%     addDigitalChannel(s,'dev5','Port0/Line0:0','OutputOnly');
    
%     %% Prepare the arduino
%     a = arduino('COM15', 'Uno');  
%     %a.writeDigitalPin("D6", 1); pause(0.005);a.writeDigitalPin("D6",0); 
%     %TODO check port, pin and pause time

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

%     % set the environment for the Time Series in PrairieView
%     tslCommand = "tsl " + envPath;
%     pl.SendScriptCommands(tslCommand);  %TODO check if this works
    
    
    %% Load Baseline variables
    % load onacid masks
    % TODO check if this loads inside of the function workspace
    load([savePath, 'redcompBMI.mat'], 'AComp');
    strcMask = obtainStrcMask(AComp, px, py);
    
    %% Create the file where to store info in case matlab crashes
    fileName = [savePath, 'bmiExp.dat'];
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    fwrite(fileID, bmiAct,'double');
    fwrite(fileID, cursor ,'double');
    fwrite(fileID, baseVector ,'double');
    fwrite(fileID, hits,'single');
    fwrite(fileID, trialStart,'single');
    fwrite(fileID, T1,'single');
    fwrite(fileID, ensemble,'single');
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',{'double',size(bmiAct),'bmiAct'; ...
        'double',size(cursor),'cursor'; ...
        'double',size(baseVector),'baseVector'; ...
        'single',size(hits),'hits'; ...
        'single',size(trialStart),'trialStart'}, 'repeat', 1); 
    m.Writable = true;
    
    %in case matlab crashes copy some info in txt
    fileID = fopen([savePath, 'bmiExp.txt'],'wt');
    fprintf(fileID,'Length %d\n',round(expectedLengthExperiment));
    fprintf(fileID,'\n%6s %12s\r\n','E1','E2');
    A = [E1; E2];
    fprintf(fileID,'%6d %12d\r\n',A);
    fclose(fileID);
    
    %initialize the values of the memmap
    m.Data.trialStart = trialStart;
    m.Data.hits = hits;

    %************************************************************************
    %*************************** RUN ********************************
    %************************************************************************

    %start the time_series scan
%     pl.SendScriptCommands("-ts");   % TODO check if this actually works 

    tic;
    while frame <= expectedLengthExperiment
        Im = pl.GetImage_2(chanIdx, px, py);
%         if Im ~= lastFrame   
%             lastFrame = Im;   % comparison and assignment takes ~4ms
            
%             % Synchronization
%             outputSingleScan(s,1);
%             pause(syncTime)
%             outputSingleScan(s,0);
            
            if nonBufferUpdateCounter == 0
                % obtain value of the neurons fluorescene
                unitVals = obtainRoi(Im, strcMask); % function to obtain Rois values
                bmiAct(:,frame) = unitVals;
                m.Data.bmiAct(:,frame) = unitVals; % 1 ms  store info
                
                % update buffer and baseval
                expHistory(:, 1: end-1) = expHistory(:, 2:end);
                expHistory(:,end) = unitVals;
                
                % calculate baseline activity and actual activity for the DFF
                signal = single(nanmean(expHistory, 2));
                if frame == initFrameBase
                    baseval = single(ones(numberNeurons,1)).*unitVals;
                elseif frame <= baseFrames
                    baseval = (baseval*(frame - 1) + signal)./frame;
                else
                    baseval = (baseval*(baseFrames - 1) + signal)./baseFrames;
                end
                baseVector(:,frame) = baseval;
                m.Data.baseVector(:,frame) = baseval; % saving in memmap
                
                % calculate of DFF
                dff = (signal - baseval) ./ baseval;
                cursor(frame) = nansum([-nansum(dff(E1)), nansum(dff(E2))]);
                disp (cursor(frame));
                m.Data.cursor(frame) = cursor(frame); % saving in memmap
                
                if BufferUpdateCounter == 0
                    % Is it a new trial?
                    if trialFlag && ~backtobaselineFlag
                        trialStart(frame) = 1;
                        m.Data.trialStart(frame) = 1;
                        trialHistory = trialHistory + 1;
                        trialFlag = 0;
                        %start running the timer again
                        disp('New Trial!')
                    end

                    if backtobaselineFlag 
                        if cursor(frame) <= back2Base 
                            backtobaselineFlag = 0;
                        end
                    else
                        if cursor(frame) >= T1      %if it hit the target
                            % VTA STIM
                            if flagVTA
                                % Arduino pulse TODO
                                a.writeDigitalPin("D6", 1); pause(syncVTA);a.writeDigitalPin("D6",0);
                                nonBufferUpdateCounter = shutterVTA;
                            end
                            % update rewardHistory
                            rewardHistory = rewardHistory + 1;
                            disp(['Trial: ', num2str(trialHistory), 'Rewards: ', num2str(rewardHistory)]);
                            % update trials and hits vector
                            hits(frame) = 1;
                            m.Data.hits(frame) = 1;
                            trialFlag = 1;
                            BufferUpdateCounter = relaxationFrames; 
                            %TODO be careful when we set the VTA about the
                            %counters
                            % relaxationtime = length of the signal buffer
                            backtobaselineFlag = 1;
                        elseif flagHolo
                            if ismember(frame, vectorHolo)
                                % holo STIM trigger goes here TODO
                            end
                        end
                    end
                else
                    BufferUpdateCounter = BufferUpdateCounter - 1;
                end
            else        
                nonBufferUpdateCounter = nonBufferUpdateCounter - 1;
            end
            frame = frame + 1;
            timeVector(frame) = toc;
%         end

    end
    pl.Disconnect();

end
% 
% % fires when main function terminates (normal, error or interruption)
function cleanMeUp(savePath, animal, day, E1, E2, T1, E1Base, E2Base)
    global pl cursor hits trialStart bmiAct baseVector timeVector btbaseVector
    disp('cleaning')
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(savePath + "BMI_online.mat", 'cursor', 'baseVector', 'hits', 'timeVector', ...
    'trialStart', 'bmiAct', 'animal', 'btbaseVector', 'day', 'E1', 'E2', 'T1', ...
    'E1Base', 'E2Base')
    pl.Disconnect();
end

