function [cursor, E2_hit, E1_hit] = ...
    dff2cursor_target_Esymm_no_constraint(dff, cal)
%NOTE: 
%dff must be of size 1xnum_neurons, or num_neuronsx1
% OUTPUT: 
% c1: cursor
% c2: E1_mean > E1_thresh
% c3: E2_subord_mean > E2_subord_thresh
%
% relevant fields from bData: 
%n_mean, n_std, decoder, E1_sel_idxs, E2_sel_idxs, E1_thresh, E2_subord_thresh, T)
%decoder: num_neurons x 1

%cond1: neural*decoder >= T
%cond2: mean E1<= E1_thresh
%cond3: mean E2 >= E2_thresh

%neural: 1 x num_neurons
%decoder: num_neurons x 1
%E_id: num_neurons x 1
%E1_thresh: 1x1
%E2_subord_thresh: 1xnum_E2
%T: 1x1

%Ensemble Identity Info:
% E1_sel = E_id==1; 
% E1_sel_idxs = find(E1_sel); 
% num_E1 = length(E1_sel_idxs); 
% 
% E2_sel = E_id==2; 
% E2_sel_idxs = find(E2_sel); 
% num_E2 = length(E2_sel_idxs); 

num_E2 = cal.neurons.num_E2; %length(bData.E2_sel_idxs); 
dff = dff(:).';
%z-score:
if(cal.cursor_zscore_bool)
    dff_z = (dff-cal.target.n_mean)./cal.target.n_std;
    dff_z = dff_z(:).'; %set dff_z to be a row
    n_analyze = dff_z;
elseif(cal.range_norm_bool)
    n_analyze = dff./cal.target.n_range; 
else
    n_analyze = dff;
end

E1 = n_analyze(cal.neurons.E1_sel_idxs); 
E1_mean = mean(E1); 
E2 = n_analyze(cal.neurons.E2_sel_idxs); 
E2_mean = mean(E2); 

cursor = n_analyze*cal.decoder;

%E2 hit: 
c1_bool = ...
    cursor >= cal.target.E2_hit_cal.T;
% c2_bool = ...
%     E2_mean >= cal.target.E2_hit_cal.E2_thresh;
% c3_bool = ...
%     E1_mean <= cal.target.E2_hit_cal.E1_thresh;

E2_hit = c1_bool; % && c2_bool && c3_bool; 

%E1 hit: 
c1_bool = ...
    cursor <= cal.target.E1_hit_cal.T;
% c2_bool = ...
%     E1_mean >= cal.target.E1_hit_cal.E1_thresh;
% c3_bool = ...
%     E2_mean <= cal.target.E1_hit_cal.E2_thresh;

E1_hit = c1_bool; % && c2_bool && c3_bool; 
