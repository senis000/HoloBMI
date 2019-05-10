function BMIAcqnvsPrairienoTrials(folder, animal, day, expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)
    %{
    Function to acquire the BMI in a prairie scope
    animal -> animal for the experiment
    day -> day for the experiment

    baselineCalibrationFile:
    %AComp_BMI: matrix for spatial filters with px*py*unit
    %zscore parameters: n_mean, n_std
    %cursor parameters: decoder
    %target parameters: 'T1', 'E1_thresh', 'E2_subord_thresh'

    vectorHolo -> vector of scheduled holo stims

    Flag description:
    flagHolosched = false;
    flagVTAsched = false;

    flagBMI: 
        determines if code detects self-generated target neural
    patterns
    flagVTAtrig: 
        determines if target neural patterns trigger VTA stim
    flagHolosched: 
        determines if holo stim will be delivered on a schedule
    flagVTAsched: 
        determines if VTA stim is delivered on a schedule

    expt_str --> Experiments:
    0) BMI
    flagBMI = true; (use self-generated hits to perform actions 
        for now actions are just to possibly send VTA stim
        in future, actions can include sending task feedback signal
    flagHolo = false; flagVTAsched = false;

    1) Pre-train (Holo -> VTA)
    flagVTAtrig = true; flagHolosched = true;
    This follows the vectorHolo schedule. If target pattern achieved,
    deliver VTA.
    flagBMI = false; flagVTAsched = false;

    2) No pre-training
    flagBMI = true; flagVTAtrig = false; flagHolosched=false;
    set experiment length to be length of pretrain + BMI
    flagVTAsched = false;
    
    3) Pre-train E3, Test E2
    baseline calibrate E3-E1, and E2-E1.  
    Pre-train E3: 
    flagVTAtrig = true; flagHolosched = true; flagBMI = false;
    flagVTAsched = false;

    Test E2:
    flagVTAtrig = true; flagHolosched = false; flagBMI = true;
    flagVTAsched = false;

    4) Pre-train orthogonal
    baseline calibrate E2-E1 shuffle, and E2-E1
    Pretrain E2-E1 shuffle:
    flagVTAtrig = true; flagHolosched = true; flagBMI = false;
    flagVTAsched = false;

    Test E2-E1: 
    flagVTAtrig = true; flagHolosched = false; flagBMI = true;
    flagVTAsched = false;

    5) Pre-train holo only
    flagVTA = false; flagHolo = true; flagBMI = false; 
    flagVTAsched = false;

    6) Random VTA
    flagVTAsched = true; flagVTAtrig = false; 
    flagBMI = false; flagHolosched = false;
    
%}
    if nargin <6
        frameRate = 30;
    end
    if nargin < 7
        vectorHolo = [];
        vectorVTA = [];
    elseif nargin == 7
        vectorVTA = [];
    end

    %%
    %**********************************************************
    %****************  PARAMETERS  ****************************
    %**********************************************************
    
    %% experiment FLAGS
    
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

    [flagBMI, flagVTAtrig, flagHolosched, flagVTAsched] = ...
        expt2bmi_flags(expt_str);
%     flagBMI       = true;
%     flagVTAtrig   = true;
%     flagHolosched = false;
%     flagVTAsched  = false;    
    
    %% BMI parameters 
    %frameRate = 30; % TODO check if it can be obtained from prairie 
    relaxationTime = 0;  % there can't be another hit in this many sec
    %back2Base = 1/2*bData.T1; % cursor must be under this value to be able to hit again

    savePath = fullfile(folder, animal, day); %[folder, animal, '/',  day, '/'];
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    % values of parameters in frames
    expectedLengthExperiment = 40*60*frameRate; % in frames
    baseFrames = round(2*60 * frameRate); % Period at the begginig without BMI to establish BL
    movingAverageFrames = 2;
    relaxationFrames = round(relaxationTime * frameRate);

    %% prairie view parameters
    chanIdx = 2; % green channel

    %% Reward/VTA parameters
    %Sound: 
    xrnd = randn(1000,1);
    reward_sound = audioplayer(xrnd, 10000); %Play sound using: play()
%     xrnd_filt = filter([1 1], 1, xrnd); 
%     reward_sound = audioplayer(xrnd_filt, 10000); %Play sound using: play()
%     play(reward_sound)
    
    %Shutter:
    flagShutter = 0; 
    if flagShutter
        shutterVTA = round(2*frameRate);
    else
        shutterVTA = 0; 
    end
    syncVTA = 0.001; % duration of the TTL
    
    %% Load BMI parameters from baseline calibration
    bData = load(fullfile(savePath, baselineCalibrationFile));
    back2Base = 1/2*bData.T1; % cursor must be under this value to be able to hit again

    %Fields: 
    %'n_mean', 'n_std',
    %'AComp_BMI', 'T1', 'decoder', 'E_id', 'E1_sel_idxs', 'E2_sel_idxs', 
    %'E1_base', 'E2_base',
    %'E2_subord_thresh', 'E1_thresh', 'E2_coeff', 'E2_subord_mean', 'E2_subord_std'

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************
    
    global pl data
%     cursor hits trialStart bmiAct baseVector timeVector %TODO remove timeVector
    
    numberNeurons = length(bData.E_id);
    
    %pre-allocating arrays
    expHistory = single(nan(numberNeurons, movingAverageFrames));  %define a windows buffer
    data.cursor = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.bmiAct = double(nan(numberNeurons, ceil(expectedLengthExperiment)));
    data.baseVector = double(nan(numberNeurons,ceil(expectedLengthExperiment)));  %define a very long vector for cursor    
    data.selfHits = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits
    data.holoHits = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.selfVTA = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.holoVTA = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.trialStart = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector trialStart
    %to debug!!! TODO REMOVE after debugging
    data.timeVector = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor

    %initializing general flags and counters 
    data.selfTargetCounter = 0; 
    data.holoTargetCounter = 0; 
    data.selfTargetVTACounter = 0; 
    data.holoTargetVTACounter = 0;
    data.schedHoloCounter = 0; 
    data.schedVTACounter = 0; 
     
    data.trialCounter = 0; %todo remove one
    trialFlag = 1;
    nonBufferUpdateCounter = 40;  %counter when we dont want to update the buffer: 
    initFrameBase = nonBufferUpdateCounter + 1;
    %beginning of experiment and VTA stim
    BufferUpdateCounter = 0;
    
    %Only useful if: flagBMI=false; flagHolosched = true; flagVTAtrig = true;
    HoloTargetWin = 10; %number of frames after a holo stim to look for target
    HoloTargetDelayTimer = 0; %if this timer is >0 check for a holo target
%     detectHoloTargetFlag = 0; %if this is 1, start looking for a holo target
    
    backtobaselineFlag = 0;
    data.frame = 1; % initialize frames
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(savePath, bData));  %in case of ctrl-c it will launch cleanmeup

%     %% Prepare the nidaq
    clear s
    s = daq.createSession('ni');
    addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
    ni_out = [0 0 0]; 
    outputSingleScan(s,ni_out);%set   
    ni_getimage = [1 0 0]; 
    ni_reward   = [0 1 0]; 
    ni_holo     = [0 0 1]; 
%       Line 0: GetImage Pulse
%       Line 1: Triggers VTA stim / reward 
%       Line 2: Triggers holo stim

    %% Prepare for Prairie
    % connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect()
    pause(2);  % pause is needed to give time to Prairie to connect

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    % Prairie commands
    pl.SendScriptCommands("-srd True 0");
    pl.SendScriptCommands("-lbs True 0");

    lastFrame = zeros(px, py); % to compare with new incoming frames

    % set the environment for the Time Series in PrairieView
    loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_40.env");
    pl.SendScriptCommands(loadCommand);   
    
    % set the path where to store the imaging data -SetSavePath (-p) "path" ["addDateTime"]
    savePrairieFiles(savePath, pl, expt_str)  
    
    
    %% load onacid masks
%     load([savePath, 'redcompBMI.mat'], 'AComp');
    strcMask = obtainStrcMask(bData.AComp_BMI, px, py);
    
    %% Create the file where to store info in case matlab crashes
    fileName = [savePath, 'bmiExp.dat'];
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    fwrite(fileID, data.cursor ,'double');
    fwrite(fileID, data.bmiAct, 'double'); 
    fwrite(fileID, data.baseVector, 'double');     
    fwrite(fileID, data.selfHits ,'single');
    fwrite(fileID, data.holoHits ,'single');
    fwrite(fileID, data.selfVTA ,'single');
    fwrite(fileID, data.holoVTA ,'single');
    fwrite(fileID, data.trialStart, 'single');    
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',...
        {'double',size(data.cursor),'cursor'; ...
        'double',size(data.bmiAct),'bmiAct'; ...
        'double',size(data.baseVector),'baseVector'; ...
        'single',size(data.selfHits),'selfHits'; ...
        'single',size(data.holoHits),'holoHits'; ...
        'single',size(data.selfVTA),'selfVTA'; ...
        'single',size(data.holoVTA),'holoVTA'; ...
        'single',size(data.trialStart),'trialStart'}, 'repeat', 1);     
    m.Writable = true;
    
%     %TODO: fix this given bData
%     %in case matlab crashes copy some info in txt
%     fileID = fopen([savePath, 'bmiExp.txt'],'wt');
%     fprintf(fileID,'Length %d\n',round(expectedLengthExperiment));
%     fprintf(fileID,'\n%6s %12s\r\n', 'bData.E1_sel_idxs', 'bData.E2_sel_idxs'); 
% %     fprintf(fileID,'\n%6s %12s\r\n', 'E1', 'E2');     
% %     A = [E1; E2];
% %     fprintf(fileID,'%6d %12d\r\n',A);
%     fclose(fileID);
    
    %initialize the values of the memmap
    m.Data.trialStart   = data.trialStart;
    m.Data.selfHits     = data.selfHits;
    m.Data.holoHits     = data.holoHits;
    m.Data.selfVTA      = data.selfVTA;
    m.Data.holoVTA      = data.holoVTA;

    %% ************************************************************************
    %*************************** RUN ********************************
    %************************************************************************

    %start the time_series scan
    pause(2); 
    pl.SendScriptCommands("-ts");  
    
    pause(2);  %empirically discovered time for the prairie to start gears
    data.frame = 1;
    tic;
    disp('STARTING RECORDING!!!')
    counterSame = 0;
    while counterSame < 500 %while data.frame <= expectedLengthExperiment
        Im = pl.GetImage_2(chanIdx, px, py);
        if ~isequal(Im,lastFrame)
            lastFrame = Im;   % comparison and assignment takes ~4ms
            outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,[0 0 0]);

            if nonBufferUpdateCounter == 0
                % obtain value of the neurons fluorescene
                unitVals = obtainRoi(Im, strcMask); % function to obtain Rois values
                data.bmiAct(:,data.frame) = unitVals;
                m.Data.bmiAct(:,data.frame) = unitVals; % 1 ms  store info
                
                % update buffer and baseval
                expHistory(:, 1: end-1) = expHistory(:, 2:end);
                expHistory(:,end) = unitVals;
                
                % calculate baseline activity and actual activity for the DFF
                signal = single(nanmean(expHistory, 2));
                if data.frame == initFrameBase
                    baseval = single(ones(numberNeurons,1)).*unitVals;
                elseif data.frame <= baseFrames
                    baseval = (baseval*(data.frame - 1) + signal)./data.frame;
                else
                    baseval = (baseval*(baseFrames - 1) + signal)./baseFrames;
                end
                data.baseVector(:,data.frame) = baseval;
                m.Data.baseVector(:,data.frame) = baseval; % saving in memmap
                
                %----------------------------------------------------------
                %Cursor
                % calculate of DFF
                dff = (signal - baseval) ./ baseval;
                [~, cursor_i, target_hit, ~] = ...
                    dff2cursor_target(dff, bData);
                data.cursor(data.frame) = cursor_i;
                m.Data.cursor(data.frame) = data.cursor(data.frame); % saving in memmap
                %----------------------------------------------------------
                
                
                if BufferUpdateCounter == 0
                    % Is it a new trial?
                    if trialFlag && ~backtobaselineFlag
                        data.trialStart(data.frame) = 1;
                        m.Data.trialStart(data.frame) = 1;
                        data.trialCounter = data.trialCounter + 1;
                        trialFlag = 0;
                        %start running the timer again
                        disp('New Trial!')
                    end

                    if backtobaselineFlag 
                        if data.cursor(data.frame) <= back2Base 
                            backtobaselineFlag = 0;
                            disp('back to baseline')
                        end
                    else
                        if target_hit      %if it hit the target
                            disp('target hit')
                            if(HoloTargetDelayTimer > 0)
                                disp('holo target was happening')
                                HoloTargetDelayTimer = 0; 
                                data.holoTargetCounter = data.holoTargetCounter + 1;
                                data.holoHits(data.frame) = 1;
                                m.Data.holoHits(data.frame) = 1;
                                
                                if flagVTAtrig
                                    disp('reward delivery!')
                                    play(reward_sound);
                                    outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out)
                                    
                                    nonBufferUpdateCounter = shutterVTA;   
                                    
                                    data.holoTargetVTACounter = data.holoTargetVTACounter+1;
                                    data.holoVTA(data.frame) = 1;
                                    m.Data.holoVTA(data.frame) = 1;
                                end
                            else
                                data.selfTargetCounter = data.selfTargetCounter + 1;
                                data.selfHits(data.frame) = 1;
                                m.Data.selfHits(data.frame) =1;
                                disp('holo target was NOT happening')
                                if(flagBMI && flagVTAtrig)
                                    nonBufferUpdateCounter = shutterVTA;
                                    disp('reward delivery! (self-target)')
                                    outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);                                    
                                    
                                    data.selfTargetVTACounter = data.selfTargetVTACounter + 1;
                                    data.selfVTA(data.frame) = 1;
                                    m.Data.selfVTA(data.frame) = 1;                                    
                                end                                
                            end

                            disp(['Trial: ', num2str(data.trialCounter), 'VTA STIMS: ', num2str(data.holoTargetVTACounter + data.selfTargetVTACounter)]);
                            % update trials and hits vector
                            trialFlag = 1;
                            BufferUpdateCounter = relaxationFrames; 
                            backtobaselineFlag = 1;
                        elseif flagHolosched
                            if ismember(data.frame, vectorHolo)
                                disp('HOLO STIM')
                                HoloTargetDelayTimer = HoloTargetWin;
                                data.schedHoloCounter = data.schedHoloCounter + 1;
                                outputSingleScan(s,ni_holo); pause(0.001); outputSingleScan(s,ni_out)
                            end
                        elseif flagVTAsched
                            if ismember(data.frame, vectorVTA)
                                disp('scheduled VTA STIM')
                                outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out)
                                nonBufferUpdateCounter = shutterVTA;                                
                                data.schedVTACounter = data.schedVTACounter + 1; 
                            end
                        end
                    end
                else
                    BufferUpdateCounter = BufferUpdateCounter - 1;
                end
            else        
                nonBufferUpdateCounter = nonBufferUpdateCounter - 1;
            end
            
            if(HoloTargetDelayTimer > 0)
                HoloTargetDelayTimer = HoloTargetDelayTimer-1;
            end
            
            data.frame = data.frame + 1;
            data.timeVector(data.frame) = toc;
            counterSame = 0;
        else
            counterSame = counterSame + 1;
        end

    end
%    pl.Disconnect();
end
% 
% % fires when main function terminates (normal, error or interruption)
function cleanMeUp(savePath, bData)
    global pl data
    disp('cleaning')
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(fullfile(savePath, "BMI_online" + datestr(datetime('now'), 'yymmddTHHMMSS') + ".mat"), 'data', 'bData')
    if pl.Connected()
        pl.Disconnect();
    end
end

