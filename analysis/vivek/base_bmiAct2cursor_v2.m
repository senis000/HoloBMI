function [est] = base_bmiAct2cursor_v2(f_preprocessed_bool, baseVector, bmiAct, movingAverageFrames, cursor_zscore_bool, bData, cursor)

%5.30.19
%Compatible with BMI code: 
%BMIAcqnvsPrairienoTrialsHoloCL_debug_enable
%done for expts from May 2019
%
%input:
% f_preprocessed_bool: if 1, f were already processed, they don't need
% to be smoothed, etc.
%baseVector: num_neurons X num_time_samples
%bmiAct: num_neurons X num_time_samples
%cursor: 1 x num_time_samples
%movingAverageFrames - how many frames to smooth for dff (non explicit
%parameter for the BMI)
%cursor_zscore_bool - whether to zscore neurons to calculate cursor (non explicit parameter for BMI)

% [dff_z, cursor, target_hit, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
%     dff2cursor_target(dff, bData, cursor_zscore_bool) 


%%
% % %Debug input: 
% f_preprocessed_bool  = 0; 
% i=1;
% baseVector      = pretrain(i).data.baseVector;
% bmiAct          = pretrain(i).data.bmiAct;
% bData           = pretrain(i).bData; 
% cursor          = pretrain(i).data.cursor;  
% % 
% % baseVector  = bmi.data.baseVector;
% % bmiAct      = bmi.data.bmiAct;
% % cursor      = bmi.data.cursor;
% % bData       = bmi.bData; 
% % 
% movingAverageFrames     = 4;
% cursor_zscore_bool      = 0 ; 

%%
%ToDo: 
%
%For pretrain, BMI, we need to convert bmiAct into neural activity for the
%cursor.  We can doublecheck by reconstructing the cursor. 
%
%if pretrain not seeded:
%number of valid samples for bmiAct and baseVector should be more than
%cursor, cuz of base updating for two minutes
%
%But if pretrain seeded, then bmiAct, baseVector, cursor should all be
%updated after "nonBufferUpdateCounter"
%%

num_neurons = length(bData.E1_base) + length(bData.E2_base); 

%%
%{
Approach:
1) preprocess the neural data
2) for the cursor valid idxs, use neural data to resimulate cursor data
%}

%%
%Preprocess the neural data: 

n_valid_idxs = find(~isnan(bmiAct(1,:))); 
est.n_valid_idxs = n_valid_idxs; 

n_valid.baseVector = ...
    baseVector(:,n_valid_idxs); 
n_valid.bmiAct = ...
    bmiAct(:,n_valid_idxs);

if(~f_preprocessed_bool)
%calculate smoothed F: 
    F_smooth = n_valid.bmiAct; 
    mean_filter = ones(movingAverageFrames, 1)/movingAverageFrames; 
    for n =1:num_neurons
        F_smooth(n,movingAverageFrames:end) = conv(n_valid.bmiAct(n,:), mean_filter, 'valid'); 
    end
    for i=1:(movingAverageFrames-1)
        F_smooth(:,i) = mean(n_valid.bmiAct(:,1:i), 2); 
    end
else
    F_smooth = bmiAct;
end
n_valid.dff = (F_smooth-n_valid.baseVector)./n_valid.baseVector;
n_valid.dff_pad = ones(size(bmiAct))*nan; 
n_valid.dff_pad(:, n_valid_idxs) = n_valid.dff; 

est.dff = n_valid.dff_pad; 


%%
%Simulate cursor on cursor valid idxs: 
c_valid_idxs = find(~isnan(bmiAct(1,:)) & ~isnan(cursor(1,:))); 
c_num_samples = length(c_valid_idxs); 
est.c_valid_idxs = c_valid_idxs; 
est.c_num_samples = c_num_samples; 

data_valid.valid_idxs = c_valid_idxs; 
data_valid.baseVector = ...
    baseVector(:,c_valid_idxs); 
data_valid.bmiAct = ...
    bmiAct(:,c_valid_idxs);
data_valid.cursor = ...
    cursor(:,c_valid_idxs); 
 
N = length(bmiAct);
est.bData           = bData; 
est.dff             = n_valid.dff_pad;
%Estimate: cursor, target_hit, c1_hit, c2_hit, c3_hit
est.data_valid      = data_valid; 
est.dff_z           = zeros(num_neurons,N); 
est.cursor          = zeros(1,N); 
est.target_hit      = zeros(1,N); 
est.c1_hit          = zeros(1,N); 
est.c2_hit          = zeros(1,N); 
est.c3_hit          = zeros(1,N); 
est.c2_val          = zeros(1,N); 
est.c3_val          = zeros(1,N); 

for i=1:c_num_samples
    a = c_valid_idxs(i); 
    [dff_z_i, cursor_i, target_hit_i, c1_bool, c2_val, c2_bool, c3_val, c3_bool] = ...
        dff2cursor_target(est.dff(:,a), bData, cursor_zscore_bool);    
    
    %Assign:
    est.dff_z(:,a)      = dff_z_i(:);
    est.cursor(a)       = cursor_i; 
    est.target_hit(a)   = target_hit_i; 
    est.c1_hit(a)          = c1_bool; 
    est.c2_hit(a)          = c2_bool; 
    est.c3_hit(a)          = c3_bool;
    est.c2_val(a)          = c2_val; 
    est.c3_val(a)          = c3_val; 
end
        
% end

% %%
% %Compare cursor estimation: 
% sum(est.cursor - data_valid.cursor > 1e-3)
% 
% %%
% h = figure;
% scatter(est.cursor, data_valid.cursor)



