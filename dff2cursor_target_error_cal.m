function [dff_z, target_hit, ...
    E2mE1, E2mE1_bool, E2mE1_error, ...
    E1_val, E1_bool, E1_error, ...
    E2_val, E2_bool, E2_error] = ...
    dff2cursor_target_error_cal(dff, cal, cursor_zscore_bool) 
%4.19.19
%INPUT: smoothed dff
%1) z-score dff
%2) calculate cursor value
%3) determine if target is hit
% OUTPUT: 
% c1: cursor
% c2: E1_mean > E1_thresh
% c3: E2_subord_mean > E2_subord_thresh
%
% relevant fields from cal: 
%n_mean, n_std, decoder, E1_sel_idxs, E2_sel_idxs, E1_thresh, E2_subord_thresh, T)
%decoder: num_neurons x 1

%E2_subord_thresh: size (num_E2-1 x 1)
%cond1: neural*decoder >= T
%cond2: E1<= E1_thresh
%cond3: E2_subord_mean >= E2_subord_thresh(E2_dom)

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

num_E2 = cal.neurons.num_E2; %length(cal.E2_sel_idxs); 

%z-score:
dff = dff(:).';
dff_z = dff; 
if(cursor_zscore_bool)
    dff_z = (dff-cal.target.n_mean)./cal.target.n_std;
    dff_z = dff_z(:).'; %set dff_z to be a row    
    n_analyze = dff_z;
else
    n_analyze = dff;
end

E1 = n_analyze(cal.neurons.E1_sel_idxs); 
E2 = n_analyze(cal.neurons.E2_sel_idxs); 

%c1: cursor
E2mE1 = n_analyze*cal.decoder; 
E2mE1_error = max(0,cal.target.T-E2mE1); 

%c2: E1_mean
E1_mean = mean(E1); 
E1_val = E1_mean;
E1_error = max(0, E1_val-cal.target.E1_thresh); 

%c3: E2_subord
E2_sum                          = sum(E2); 
[E2_dom_samples, E2_dom_sel]    = max(E2, [], 2); %E2_dom  CAREFUL THIS MAY BRING TWO!!
% E2_dom_samples = n_analyze(E2_dom_sel(1));
E2_subord_mean                  = (E2_sum - E2_dom_samples)/(num_E2-1);
E2_val                          = E2_subord_mean;
E2_error                        = max(0, cal.target.E2_subord_thresh(E2_dom_sel)-E2_val); 

%Booleans:
E2mE1_bool = E2mE1 >= cal.target.T; 
E1_bool = E1_mean <= cal.target.E1_thresh;
E2_bool = E2_subord_mean >= cal.target.E2_subord_thresh(E2_dom_sel);

%target hit
target_hit = E2mE1_bool & E1_bool & E2_bool; 



% %1) (E2-E1) >= T
% % which means T-(E2-E1) <= 0
% E2mE1_error     = max(0, cursor_target - cursor_obs); 
% 
% %2) E1 <= mu
% % which means E1-mu <= 0
% E1_error        = max(0, E1_mean_analyze - E1_thresh); 
% 
% %3) E2_subord > mu
% % which means mu-E2_subord 
% 
% E2_error        = max(0, E2_subord_thresh(E2_dom_sel) - E2_subord_mean_analyze); 