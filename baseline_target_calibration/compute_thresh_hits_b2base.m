function [E_hit_cal, E_hit_data] = compute_thresh_hits_b2base(E1_bool, E_cal, ...
    cursor_obs, ...
    mean_E2, mean_E1, ...
    b2base_coeff, b2baseFrameThresh)
%NOTE:
%This code doesn't use 'mean_E2', 'mean_E1' to calculate the hits.  Can
%input []
%
%Return Values:
% hit_cal
% hit_cal.T                   = T_value;
% hit_cal.num_E2_valid        = num_E2_valid; 
% hit_cal.num_E1_valid        = num_E1_valid; 
% hit_cal.num_hits_b2base     = num_hits_b2base; 
% hit_cal.num_hits_no_b2base  = num_hits_no_b2base; 
% hit_cal.reward_prob_per_frame   = reward_prob_per_frame; 
%
% hit_data
% hit_data                    = hit_cal;
% hit_data.T_hits             = T_hits; 
% hit_data.E2_valid           = E2_valid; 
% hit_data.E1_valid           = E1_valid; 
% hit_data.T_idxs_no_b2base   = T_idxs_no_b2base;
% hit_data.T_idxs_b2base      = T_idxs_b2base;

T_value = E_cal.T;
b2base_thresh                   = b2base_coeff*T_value;

if(E1_bool)
%     T_neg                           = prctile(-cursor_obs, T_prctile); 
%     T_value                         = -T_neg; 
%     b2base_thresh                   = b2base_coeff*T_value;        
    T_hits                          = find(cursor_obs <= T_value); 
    E2_valid                        = find(mean_E2 <= E_cal.E2_thresh); 
    E1_valid                        = find(mean_E1 >= E_cal.E1_thresh);    
else
    T_hits                          = find(cursor_obs >= T_value); 
    E2_valid                        = find(mean_E2 >= E_cal.E2_thresh); 
    E1_valid                        = find(mean_E1 <= E_cal.E1_thresh);
end

%Intersection of these are the target hits without b2base constraint:
T_idxs_no_b2base        = T_hits; %intersect(intersect(E1_valid, E2_valid), T_hits); 
 
hits_valid              = ones(length(T_idxs_no_b2base),1); 
if length(T_idxs_no_b2base) > 1
    for i = 2:length(T_idxs_no_b2base)
        if(E1_bool)
            b2base_bool = sum(cursor_obs(T_idxs_no_b2base(i-1):T_idxs_no_b2base(i)) >= b2base_thresh) >= b2baseFrameThresh; 
        else
            b2base_bool = sum(cursor_obs(T_idxs_no_b2base(i-1):T_idxs_no_b2base(i)) <= b2base_thresh) >= b2baseFrameThresh; 
        end
        hits_valid(i) = b2base_bool; 
    end
end
T_idxs_b2base               = T_idxs_no_b2base(find(hits_valid));    

num_E2_valid                = length(E2_valid); 
num_E1_valid                = length(E1_valid); 
num_hits_b2base             = length(T_idxs_b2base); 
num_hits_no_b2base          = length(T_idxs_no_b2base); 
reward_prob_per_frame       = ...
    num_hits_b2base/length(cursor_obs);

%ASSIGN:
%hit_cal
E_hit_cal                       = E_cal; 
E_hit_cal.T                   = T_value;
E_hit_cal.b2base_thresh       = b2base_thresh; 
E_hit_cal.num_E2_valid        = num_E2_valid; 
E_hit_cal.num_E1_valid        = num_E1_valid; 
E_hit_cal.num_hits_b2base     = num_hits_b2base; 
E_hit_cal.num_hits_no_b2base  = num_hits_no_b2base; 
E_hit_cal.reward_prob_per_frame   = reward_prob_per_frame; 


%hit_data
E_hit_data                    = E_hit_cal;
E_hit_data.T_hits             = T_hits; 
E_hit_data.E2_valid           = E2_valid; 
E_hit_data.E1_valid           = E1_valid; 
E_hit_data.T_idxs_no_b2base   = T_idxs_no_b2base;
E_hit_data.T_idxs_b2base      = T_idxs_b2base;    

end