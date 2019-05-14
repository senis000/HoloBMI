function [cursor_vec, hit_vec, num_hits, c1_vec, c2_vec, c3_vec] = ...
    bmi_sim(f, frameRate, bData, cursor_zscore_bool)
%Input: 
% dff: num_samples X num_neurons

numberNeurons = size(f,2); 
num_samples = size(f,1); 

Fbuffer = single(nan(numberNeurons, movingAverageFrames));  %define a windows buffer
data.trialCounter = 0; %todo remove one
trialFlag = 1;

baseFrames = 2*60*frameRate; 
nonBufferUpdateCounter = 40;  %counter when we dont want to update the buffer: 
initFrameBase = nonBufferUpdateCounter + 1;
relaxationFrames = 0; %number of frames which must pass 
%beginning of experiment and VTA stim
BufferUpdateCounter = 0;
backtobaselineFlag = 0;

%OUTPUTS: 
c1_vec = zeros(num_samples,1)+nan; 
c2_vec = zeros(num_samples,1)+nan;
c3_vec = zeros(num_samples,1)+nan;
hit_vec = zeros(num_samples,1); 
num_hits = 0; 
cursor_vec = zeros(num_samples,1)+nan; 

baseBuffer_full = 0; %bool indicating the Fbuffer filled
for frame_i = 1:num_samples
    if nonBufferUpdateCounter == 0
        % obtain value of the neurons fluorescence
        unitVals = f(frame_i,:).'; 

        % update F buffer
        Fbuffer(:, 1:end-1) = Fbuffer(:, 2:end);
        Fbuffer(:,end) = unitVals;

        % calculate F0 baseline activity 
        if frame_i == initFrameBase
    %                     baseval = single(ones(numberNeurons,1)).*unitVals;
            baseval = single(ones(numberNeurons,1)).*unitVals/baseFrames;
            %---
        elseif frame_i <= (initFrameBase+baseFrames)
    %                     baseval = base(baseval*(frame_i - 1) + signal)./frame_i;
            baseval = baseval + unitVals/baseFrames;
            disp(frame_i);
            if frame_i == (initFrameBase+baseFrames)
                baseBuffer_full = 1;
%                 disp('baseBuffer FULL!'); 
            end
        elseif frame_i > (initFrameBase+baseFrames)
            baseval = (baseval*(baseFrames - 1) + unitVals)./baseFrames;
        end

        %Smooth F
        Fsmooth = single(nanmean(Fbuffer, 2));                

        if baseBuffer_full
            %----------------------------------------------------------
            %Cursor
            % calculate (smoothed) DFF
            dff = (Fsmooth - baseval) ./ baseval;
            %Passing smoothed dff to "decoder"
            [~, cursor_i, target_hit, c1_bool, ~, c2_bool, ~, c3_bool] = ...
                dff2cursor_target(dff, bData, cursor_zscore_bool);
    %                     data.bmidffz(:,frame_i) = dff_z;
            cursor_vec(frame_i) = cursor_i; 
            c1_vec(frame_i) = c1_bool; 
            c2_vec(frame_i) = c2_bool;
            c3_vec(frame_i) = c3_bool;
        end
        
        if BufferUpdateCounter == 0 && baseBuffer_full
            % Is it a new trial?
            if trialFlag && ~backtobaselineFlag
                data.trialCounter = data.trialCounter + 1;
                trialFlag = 0;
%                 disp('New Trial!');
            end        
            if backtobaselineFlag 
                if data.cursor(data.frame) <= back2Base 
                    backtobaselineFlag = 0;
%                     disp('back to baseline')
                end
            else
                if target_hit      %if it hit the target
%                     disp('target hit')
                    %Self hit:
                    hit_vec(frame_i) = 1;
                    num_hits = num_hits+1;                                   

                    BufferUpdateCounter = relaxationFrames; 
                    backtobaselineFlag = 1;
%                     disp(['Trial: ', num2str(data.trialCounter), 'VTA STIMS: ', num2str(data.holoTargetVTACounter + data.selfTargetVTACounter)]);
                    % update trials and hits vector
                    trialFlag = 1;                                    
                end
            end
        else
            BufferUpdateCounter = BufferUpdateCounter - 1;
        end
    else        
        nonBufferUpdateCounter = nonBufferUpdateCounter - 1;
    end
end