function [saveFile] = BMI_CLDA(folder, animal, day, ...
    expt_str, cal_init, task_settings, a, vectorHolo, vectorVTA, ...
    base_cursor_samples, ...
    debug_bool, debug_input)
% BMIAcqnvsPrairienoTrialsHoloCL_fb_debug_enable(folder, animal, day, ...
%     expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
%     cursor_zscore_bool, debug_bool, debug_input)
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
    
    global cal
    cal = cal_init; 
    %% BMI parameters 
    savePath = fullfile(folder, animal, day); %[folder, animal, '/',  day, '/'];
    saveFile = fullfile(savePath, ['CLDA_', datestr(datetime('now'), 'yymmddTHHMMSS'), '.mat']);
    
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    reward_enable           = 0; 
    frameRate               = task_settings.frameRate; % TODO check if it can be obtained from prairie 
    relaxationTime          = task_settings.relaxationTime;  % there can't be another hit in this many sec
    dilation_factor         = 1.2; 
    expectedLengthExperiment = ceil(task_settings.bmi_len*frameRate*dilation_factor); % in frames
    %EDIT HERE
    prefixFrames            = task_settings.prefix_win; 
    baseFrames              = task_settings.f0_win; 
    % Period at the beginning to establish f0 baseline without BMI
    movingAverageFrames     = task_settings.dff_win;
    relaxationFrames        = round(relaxationTime * frameRate);
    E1_back2Base            = cal.target.E1_hit_cal.b2base_thresh;
    E2_back2Base            = cal.target.E2_hit_cal.b2base_thresh;
%     back2Base               = cal.target.E2_hit_cal.T*task_settings.b2base_coeff;
    %In order to hit target again, cursor must be under this value for at least 
    rewardDelayFrames       = task_settings.rewardDelayFrames;
    %Number of frames between target achievement and reward delivery
    back2BaseFrameThresh    = task_settings.back2BaseFrameThresh; 
    %need to be back2Base for two frames before another target can be achieved
    
    %% CLDA Parameters
    num_cursor_base         = length(base_cursor_samples); 
    num_cursor_online       = 0; %Needs to be updated
    num_cursor_valid        = num_cursor_base + num_cursor_online; 
    last_clda               = 0; %Needs to be updated with 'num_cursor_online'
    %Perform CLDA (Update parameters) when: 
    %num_cursor_valid >= task_settings.clda.min_samples_for_adapt
    %&&
    %num_cursor_online - last_clda >=
    %task_settings.clda.period_for_adapt
   
    
    
    %% prairie view parameters
    chanIdx = 2; % green channel

    %% Reward/VTA parameters
    %Sound: 
%     xrnd = randn(1000,1);
%     reward_sound = audioplayer(xrnd, 10000); %Play sound using: play()
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

    %%
    %*********************************************************************
    %******************  INITIALIZE  ***********************************
    %*********************************************************************
    
    global pl data
%     cursor hits trialStart bmiAct baseVector timeVector %TODO remove timeVector
    
    numberNeurons   = cal.neurons.num_neurons;%length(bData.E_id);
    
    %pre-allocating arrays
    Fbuffer         = single(nan(numberNeurons, movingAverageFrames));  %define a windows buffer
    data.cursor     = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.fb_freq    = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for fb_freq
    data.bmiAct     = double(nan(numberNeurons, ceil(expectedLengthExperiment)));
%     data.bmidffz = double(nan(numberNeurons, ceil(expectedLengthExperiment)));
    data.baseVector = double(nan(numberNeurons,ceil(expectedLengthExperiment)));  %define a very long vector for cursor    
    data.selfHits   = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits
    data.E1Hits     = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits
    data.holoHits   = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.selfVTA    = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.holoVTA    = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.holoDelivery   = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector for hits    
    data.trialStart = single(zeros(1,ceil(expectedLengthExperiment)));  %define a very long vector trialStart
    %to debug!!! TODO REMOVE after debugging
    data.timeVector = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.vectorHolo = vectorHolo; 
    data.vectorHoloCL   = vectorHolo;     
    data.vectorVTA  = vectorVTA; 
    %CLDA: 
    data.E2_T       = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.E1_T       = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    data.mid_T      = double(nan(1,ceil(expectedLengthExperiment)));  %define a very long vector for cursor
    
    %initializing general flags and counters 
    data.E1TargetCounter = 0;
    data.selfTargetCounter = 0; 
    data.holoTargetCounter = 0; 
    data.selfTargetVTACounter = 0; 
    data.holoTargetVTACounter = 0;
    data.schedHoloCounter = 0; 
    data.schedVTACounter = 0; 
     
    data.trialCounter = 0; %todo remove one
    trialFlag = 1;
    nonBufferUpdateCounter = prefixFrames;  %counter when we dont want to update the buffer
    %initialize nonBufferUpdateCounter to 'prefixFrames' in order to
    %exclude these frames when BMI is starting
    initFrameBase = nonBufferUpdateCounter + 1;
    %beginning of experiment and VTA stim
    BufferUpdateCounter = 0;
    
    %Only useful if: flagBMI=false; flagHolosched = true; flagVTAtrig = true;
    HoloTargetWin = 20; %number of frames after a holo stim to look for target
    HoloTargetDelayTimer = 0; %if this timer is >0 check for a holo target
%     detectHoloTargetFlag = 0; %if this is 1, start looking for a holo target

    deliver_reward      = 0; 
    rewardDelayCounter  = 0;  
    
    back2BaseCounter = 0;
    %Counts how many frames the cursor is in baseline range, if
    %'backtobaselineFlag' is on
    %TODO: make this an input/variable loaded from calibration
    E1_backtobaselineFlag = 0;
    E2_backtobaselineFlag = 0;
    
    data.frame = 1; % initialize frames
    
    %% Cleaning 
    finishup = onCleanup(@() cleanMeUp(saveFile, cal_init, task_settings, debug_bool));  %in case of ctrl-c it will launch cleanmeup

%     %% Prepare the nidaq
    if(~debug_bool)
        clear s
        s = daq.createSession('ni');
        addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
        ni_out = [0 0 0]; 
        outputSingleScan(s,ni_out);%set   
        ni_getimage = [1 0 0]; 
        ni_reward   = [0 1 0]; 
        ni_holo     = [0 0 1]; 
    end
%       Line 0: GetImage Pulse
%       Line 1: Triggers VTA stim / reward 
%       Line 2: Triggers holo stim
%TODO: line for audio feedback ts

    %% Prepare for Prairie
    % connection to Prairie
    if(~debug_bool)
        pl = actxserver('PrairieLink.Application');
        pl.Connect()
        pause(2);  % pause is needed to give time to Prairie to connect

        % Prairie variables
        px = pl.PixelsPerLine();
        py = pl.LinesPerFrame();

        % Prairie commands
        pl.SendScriptCommands('-srd True 0');
        pl.SendScriptCommands('-lbs True 0');

        lastFrame = zeros(px, py); % to compare with new incoming frames

        % set the environment for the Time Series in PrairieView
        %TODO: don't hard code this, take it from the settings: 
        loadCommand = ['-tsl ' task_settings.clda_env]
        pl.SendScriptCommands(loadCommand);   

        % set the path where to store the imaging data -SetSavePath (-p) "path" ["addDateTime"]
        savePrairieFiles(savePath, pl, expt_str)  
    else
        px = 512; 
        py = 512; 
        lastFrame = zeros(px, py); % to compare with new incoming frames
        Im = zeros(px, py); 
    end
    
    
    %% load ROI masks
    %Data was saved in 'redcompBMI.mat' in 'baseline2target_vE1strict.m'
    %TODO: just pass the strcMask
    if(~debug_bool)
%         load(fullfile(savePath, 'redcompBMI.mat'), 'strcMask');
        load(cal.paths.BMI_roi_path, 'strcMask');
    end
    
    %% Create the file where to store info in case matlab crashes
    fileName = [savePath, 'bmiExp.dat'];
    % creates a file with the correct shape
    fileID = fopen(fileName,'w');
    fwrite(fileID, data.cursor ,'double');
    fwrite(fileID, data.fb_freq,'double');
    fwrite(fileID, data.bmiAct, 'double'); 
    fwrite(fileID, data.baseVector, 'double');     
    fwrite(fileID, data.selfHits ,'single');
    fwrite(fileID, data.E1Hits ,'single');
    fwrite(fileID, data.holoHits ,'single');
    fwrite(fileID, data.selfVTA ,'single');
    fwrite(fileID, data.holoVTA ,'single');
    fwrite(fileID, data.trialStart, 'single'); 
    fwrite(fileID, data.E2_T, 'double'); 
    fwrite(fileID, data.E1_T, 'double'); 
    fwrite(fileID, data.mid_T, 'double'); 
    
    fclose(fileID);
    % maps the file into memory
    m = memmapfile(fileName, 'Format',...
        {'double',size(data.cursor),'cursor'; ...
        'double',size(data.fb_freq),'fb_freq'; ...
        'double',size(data.bmiAct),'bmiAct'; ...
        'double',size(data.baseVector),'baseVector'; ...
        'single',size(data.selfHits),'selfHits'; ...
        'single',size(data.E1Hits),'E1Hits'; ...
        'single',size(data.holoHits),'holoHits'; ...
        'single',size(data.selfVTA),'selfVTA'; ...
        'single',size(data.holoVTA),'holoVTA'; ...
        'single',size(data.trialStart),'trialStart'; ...
        'double',size(data.E2_T), 'E2_T';...
        'double',size(data.E1_T), 'E1_T'; ...
        'double',size(data.mid_T), 'mid_T'}, 'repeat', 1);
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
    m.Data.E1Hits       = data.E1Hits;
    m.Data.holoHits     = data.holoHits;
    m.Data.selfVTA      = data.selfVTA;
    m.Data.holoVTA      = data.holoVTA;

    %%
    
    %% ************************************************************************
    %*************************** RUN ********************************
    %************************************************************************

    %start the time_series scan
    if(~debug_bool)
        pause(2); 
        pl.SendScriptCommands('-ts');  
        pause(2);  %empirically discovered time for the prairie to start gears
    end
    data.frame = 1;
    
    disp('STARTING RECORDING!!!')
    counterSame = 0; %Counts how many frames are the same as past,
    counterSameThresh = 500; %TODO put in task_settings
    baseBuffer_full = 0; %bool indicating the Fbuffer filled
    %---
    disp('baseBuffer filling!...')
    while (~debug_bool && counterSame < counterSameThresh) || (debug_bool && data.frame < size(debug_input,2)) %while data.frame <= expectedLengthExperiment
        if ~debug_bool
            Im = pl.GetImage_2(chanIdx, px, py);
        else
            Im = zeros(px, py); 
        end
        
        if ~isequal(Im,lastFrame) || debug_bool
            tic; %Start timing to see length of an iteration
            if(~debug_bool)
                lastFrame = Im;   % comparison and assignment takes ~4ms
                outputSingleScan(s,ni_getimage); pause(0.001); outputSingleScan(s,[0 0 0]);
            end
            
            if nonBufferUpdateCounter == 0
                % obtain value of the neurons fluorescence
                if(~debug_bool)
                    unitVals = obtainRoi(Im, strcMask); % function to obtain Rois values
                else
                    unitVals = debug_input(:,data.frame); 
                end
                data.bmiAct(:,data.frame) = unitVals;
                m.Data.bmiAct(:,data.frame) = unitVals; % 1 ms  store info
                
                % update F buffer
                Fbuffer(:, 1:end-1) = Fbuffer(:, 2:end);
                Fbuffer(:,end) = unitVals;
                
                % calculate F0 baseline activity 
                if data.frame == initFrameBase
%                     baseval = single(ones(numberNeurons,1)).*unitVals;
                    baseval = single(ones(numberNeurons,1)).*unitVals/baseFrames;
                    %---
                elseif data.frame <= (initFrameBase+baseFrames)
%                     baseval = base(baseval*(data.frame - 1) + signal)./data.frame;
                    baseval = baseval + unitVals/baseFrames;
                    disp(data.frame);
                    if data.frame == (initFrameBase+baseFrames)
                        baseBuffer_full = 1;
                        disp('baseBuffer FULL!'); 
                    end
                elseif data.frame > (initFrameBase+baseFrames)
                    baseval = (baseval*(baseFrames - 1) + unitVals)./baseFrames;
                end
                data.baseVector(:,data.frame) = baseval;
                m.Data.baseVector(:,data.frame) = baseval; % saving in memmap
                
                %Smooth F
                Fsmooth = single(nanmean(Fbuffer, 2));                
                
                if baseBuffer_full
                    %----------------------------------------------------------
                    %Cursor
                    % calculate (smoothed) DFF
                    dff = (Fsmooth - baseval) ./ baseval;
                    %Passing smoothed dff to "decoder"
%                     [cursor_i, target_hit, c1_bool, ~, c2_bool, ~, c3_bool] = ...
%                         dff2cursor_target_v2(dff, cal);
                    [cursor_i, E2_hit, E1_hit] = ...
                        dff2cursor_target_Esymm_no_constraint(dff, cal);
%                     [cursor_i, E2_hit, E1_hit] = ...
%                         dff2cursor_target_Esymm(dff, cal); 
                    target_hit = E2_hit; 
%                     data.bmidffz(:,data.frame) = dff_z;
%--------------------------------------------------------------------------
                    data.cursor(data.frame) = cursor_i;
                    m.Data.cursor(data.frame) = data.cursor(data.frame); % saving in memmap
                    num_cursor_online   = num_cursor_online + 1; 
                    num_cursor_valid    = num_cursor_base + num_cursor_online; 
%--------------------------------------------------------------------------
                    disp(['Cursor: ' num2str(cursor_i)]); 

                    %fb: 
%--------------------------------------------------------------------------
                    fb_freq_i = cursor2audio_freq_middle_match(cursor_i, cal);
%                     fb_freq_i = cursor2audio_freq(cursor_i, cal);  
%                     disp(['FB Freq: ' num2str(fb_freq_i)]);
                    data.fb_freq(data.frame) = fb_freq_i;
                    m.Data.fb_freq(data.frame) = data.fb_freq(data.frame); % saving in memmap    
%--------------------------------------------------------------------------                    
                    if(task_settings.fb.fb_bool && ~debug_bool)
                        %Send tone arduino
                        playTone(a,...
                            task_settings.fb.arduino.pin,...
                            fb_freq_i,...
                            task_settings.fb.arduino.duration)
                    end
                    
%                     disp(['Target : ' num2str(target_hit)]); 
%                     disp(['C1 - cursor: ' num2str(c1_bool)]); 
%                     disp(['C2 - E1 : ' num2str(c2_bool)]); 
%                     disp(['C3 - E2 subord : ' num2str(c3_bool)]);                     
                    % c1: cursor
                    % c2: E1_mean > E1_thresh
                    % c3: E2_subord_mean > E2_subord_thresh                    
                    %----------------------------------------------------------
                end
                
                if (BufferUpdateCounter == 0) && baseBuffer_full
%                     disp('HERE'); 
                    % Is it a new trial?
                    if trialFlag && ~(E1_backtobaselineFlag || E2_backtobaselineFlag)
                        data.trialStart(data.frame) = 1;
                        m.Data.trialStart(data.frame) = 1;
                        data.trialCounter = data.trialCounter + 1;
                        trialFlag = 0;
                        %start running the timer again
                        disp('New Trial!')
                    end 
                    if (E1_backtobaselineFlag || E2_backtobaselineFlag)
                        if(E1_backtobaselineFlag)
                            if data.cursor(data.frame) >= E1_back2Base 
                                back2BaseCounter = back2BaseCounter+1;
                            end
                        elseif(E2_backtobaselineFlag)
                            if data.cursor(data.frame) <= E2_back2Base 
                                back2BaseCounter = back2BaseCounter+1;
                            end
                        end
                        
                        if back2BaseCounter >= back2BaseFrameThresh
                            E1_backtobaselineFlag = 0; 
                            E2_backtobaselineFlag = 0; 
                            back2BaseCounter = 0;
                            disp('back to baseline')
                        end
                    else
%                         disp('HERE2'); 
                        if E1_hit
                            %Update count
                            data.E1TargetCounter = data.E1TargetCounter + 1;
                            %Save the point it happened at
                            data.E1Hits(data.frame) = 1;
                            m.Data.E1Hits(data.frame) = 1;
                            disp('E1 Hit!'); 
                            disp(['Num E1 Hits: ', num2str(data.E1TargetCounter)]); 
                            E1_backtobaselineFlag = 1; 
                        elseif target_hit      %if it hit the target
                            disp('target hit')
                            
                            %Holo Triggered: 
                            %----------------------------------------------
                            if(HoloTargetDelayTimer > 0)
                                disp('Holo Target Achieved')
                                HoloTargetDelayTimer = 0; 
                                data.holoTargetCounter = data.holoTargetCounter + 1;
                                data.holoHits(data.frame) = 1;
                                m.Data.holoHits(data.frame) = 1;
                                
                                if flagVTAtrig
                                    disp('RewardTone delivery!')
                                    if(~debug_bool)
                                        play(reward_sound);
%                                         outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out)
                                    end
                                    rewardDelayCounter = rewardDelayFrames; 
                                    deliver_reward = 1;                                     
                                    nonBufferUpdateCounter = shutterVTA;   
                                    
                                    data.holoTargetVTACounter = data.holoTargetVTACounter+1;
                                    data.holoVTA(data.frame) = 1;
                                    m.Data.holoVTA(data.frame) = 1;
                                end
                                
                                %Back to baseline, and new trial
                                BufferUpdateCounter = relaxationFrames; 
                                E2_backtobaselineFlag = 1;
                                disp(['Trial: ', num2str(data.trialCounter), 'VTA STIMS: ', num2str(data.holoTargetVTACounter + data.selfTargetVTACounter)]);
                                % update trials and hits vector
                                trialFlag = 1;
                                
                            %Self Hit!
                            %----------------------------------------------                                
                            else
                                %Self hit:
                                data.selfTargetCounter = data.selfTargetCounter + 1;
                                data.selfHits(data.frame) = 1;
                                m.Data.selfHits(data.frame) =1;
                                disp('self hit')
                                if(flagBMI && flagVTAtrig)
                                    nonBufferUpdateCounter = shutterVTA;
                                    disp('Target Achieved! (self-target)')
%                                     disp('RewardTone delivery!')
%                                     if ~debug_bool
%                                         play(reward_sound);
%                                     end                                        
                                    rewardDelayCounter = rewardDelayFrames; 
                                    deliver_reward = 1;                                         
%                                         outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);                                    

                                    
                                    data.selfTargetVTACounter = data.selfTargetVTACounter + 1;
                                    data.selfVTA(data.frame) = 1;
                                    m.Data.selfVTA(data.frame) = 1;                                    

                                    BufferUpdateCounter = relaxationFrames; 
                                    E2_backtobaselineFlag = 1;
                                    disp(['Trial: ', num2str(data.trialCounter), 'VTA STIMS: ', num2str(data.holoTargetVTACounter + data.selfTargetVTACounter)]);
                                    % update trials and hits vector
                                    trialFlag = 1;                                    
                                else
                                    disp(['Num Self Hits: ', num2str(data.selfTargetCounter)]); 
                                end
                            end
                        end
                        if ~trialFlag
%                             disp(['HERE ' num2str(data.frame)]); 

                            %Scheduled Stimulation
                            %----------------------------------------------
                            if flagHolosched
                                if ismember(data.frame, data.vectorHoloCL)
                                    disp('SCHEDULED HOLO STIM'); 
                                    currHoloIdx = find(data.vectorHoloCL == data.frame);                                                                   
                                        %Check E1, if lower than threshold, do
                                        %stim, and save frame
                                    if(c2_bool)                                    
                                        disp('HOLO STIM')
                                        HoloTargetDelayTimer = HoloTargetWin;
                                        data.schedHoloCounter = data.schedHoloCounter + 1;
                                        if(~debug_bool)
                                            outputSingleScan(s,ni_holo); pause(0.001); outputSingleScan(s,ni_out)                                    
                                        end
                                        %Also, save the frame we do this!!
                                        data.holoDelivery(data.frame) = 1;
                                    else
                                        data.vectorHoloCL(currHoloIdx:end) = data.vectorHoloCL(currHoloIdx:end)+1;
                                    end
                                end
                            elseif flagVTAsched
                                if ismember(data.frame, vectorVTA)
                                    disp('scheduled VTA STIM')
                                    disp('RewardTone delivered!'); 
                                    if(~debug_bool)
                                        play(reward_sound); 
%                                         outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out)
                                    end
                                    rewardDelayCounter = rewardDelayFrames; 
                                    deliver_reward = 1; 

                                    nonBufferUpdateCounter = shutterVTA;                                
                                    data.schedVTACounter = data.schedVTACounter + 1; 
                                end
                            end
                        end
                    end
                else
                    if(BufferUpdateCounter>0)
                        BufferUpdateCounter = BufferUpdateCounter - 1;
                    end
                end
            else
                if(nonBufferUpdateCounter>0)
                    nonBufferUpdateCounter = nonBufferUpdateCounter - 1;
                end
            end
            
            if(HoloTargetDelayTimer > 0)
                HoloTargetDelayTimer = HoloTargetDelayTimer-1;
            end
            
            if(rewardDelayCounter > 0)
                rewardDelayCounter = rewardDelayCounter -1; 
            elseif(reward_enable && deliver_reward && rewardDelayCounter==0)
                if(~debug_bool)
                    outputSingleScan(s,ni_reward); pause(0.001); outputSingleScan(s,ni_out);
                    %This triggers arduino to control solenoid
                end
                deliver_reward = 0; 
                disp('reward delivered!'); 
            end
            
            %Check if CLDA should be done: 
            if(num_cursor_valid > task_settings.clda.min_samples_for_adapt && ...
                    (num_cursor_online - last_clda) > task_settings.clda.period_for_adapt)
                disp('Adapting Thresholds!'); 
                last_clda = num_cursor_online; 
                adapt_params(task_settings, base_cursor_samples, data.frame);
                data.E1_T(data.frame)   = cal.target.E1_hit_cal.T; 
                data.E2_T(data.frame)   = cal.target.E2_hit_cal.T; 
                data.mid_T(data.frame)  = cal.target.cursor_middle; 
                disp('E1 T: ')
                cal.target.E1_hit_cal.T
                disp('E2 T: ')
                cal.target.E2_hit_cal.T                
                
                E1_back2Base = cal.target.E1_hit_cal.b2base_thresh;
                E2_back2Base = cal.target.E2_hit_cal.b2base_thresh;
            end
            
            %Updates:
            data.frame = data.frame + 1;
            data.timeVector(data.frame) = toc;
            counterSame = 0;
            if data.timeVector(data.frame) < 1/(frameRate*1.2)
                pause(1/(frameRate*1.2) - data.timeVector(data.frame))
            end
            
        else
            counterSame = counterSame + 1;
        end
    end
%    pl.Disconnect();
%     save(fullfile(savePath, ['BMI_online', datestr(datetime('now'), 'yymmddTHHMMSS'), '.mat']), 'data', 'bData')
end
% 
% % fires when main function terminates (normal, error or interruption)
function cleanMeUp(saveFile, cal_init, task_settings, debug_bool)
    global pl data cal
    disp('cleaning')
    % evalin('base','save baseVars.mat'); %do we want to save workspace?
    % saving the global variables
    save(saveFile, 'data', 'cal', 'task_settings', 'cal_init')
    if ~debug_bool
        if pl.Connected()
            pl.Disconnect();
        end
    end
end

function adapt_params(task_settings, base_cursor_samples, frame)
    global data cal
    
    %Get the valid online cursor samples
    first_valid_idx = 1 +  ...
        task_settings.prefix_win + ...
        task_settings.f0_win; 
    cursor_valid = data.cursor(first_valid_idx:frame); 
    
    %Join with the base_cursor_samples: 
    cursor_samples = [base_cursor_samples(:); cursor_valid(:)]; 
    %Choose the samples: 
    if task_settings.clda.use_win_for_adapt
        first_idx_use = length(cursor_samples)-task_settings.clda.win_for_adapt+1; 
        cursor_samples_use = cursor_samples(first_idx_use:end); 
    else
        cursor_samples_use = cursor_samples; 
    end
    
    %Update the percentiles: 
    E1_T    = prctile(cursor_samples_use, 100-cal.target.T_prctile); 
    E2_T    = prctile(cursor_samples_use, cal.target.T_prctile); 
    mid_T   = prctile(cursor_samples_use, task_settings.fb.middle_prctile); 
    
%     cal_CLDA                                    = cal;
    cal.target.E1_hit_cal.T                = E1_T;
    cal.target.E1_hit_cal.b2base_thresh    = E1_T*task_settings.b2base_coeff;
    cal.target.E2_hit_cal.T                = E2_T; 
    cal.target.E2_hit_cal.b2base_thresh    = E2_T*task_settings.b2base_coeff;
    cal.target.cursor_middle               = mid_T; 
    
    %Update the feedback parameters: 
    [cal] = cal2fb_Tsymm_midmapped(cal, task_settings);
end

