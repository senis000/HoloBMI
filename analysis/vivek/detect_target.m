function [target_hit, c1_val, c1_bool, c2_val, c2_bool, c3_val, c3_bool, c3_thresh] = ...
    detect_target(neural, decoder, E_id, E1_thresh, E2_subord_thresh, T)
%4.17.19
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
E1_sel = E_id==1; 
E1_sel_idxs = find(E1_sel); 
num_E1 = length(E1_sel_idxs); 

E2_sel = E_id==2; 
E2_sel_idxs = find(E2_sel); 
num_E2 = length(E2_sel_idxs); 
 
E1 = neural(E1_sel_idxs); 
E2 = neural(E2_sel_idxs); 


%c1: cursor
cursor = neural*decoder;
c1_val = cursor; 

%c2: E1_mean
E1_mean = mean(E1); 
c2_val = E1_mean;

%c3: E2_subord
E2_sum = sum(E2); 
[E2_dom, E2_dom_sel]    = max(E2, [], 2);
E2_subord_mean          = (E2_sum - E2_dom_samples)/(num_E2-1);
c3_val = E2_subord_mean;

c1_bool = cursor >= T; 
c2_bool = E1_mean <= E1_thresh;
c3_bool = E2_subord_mean >= E2_subord_thresh(E2_dom_sel);

%arget hit>
target_hit = c1_bool & c2_bool & c3_bool; 