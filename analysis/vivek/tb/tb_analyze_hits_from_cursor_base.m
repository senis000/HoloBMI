%tb_analyze_hits_from_cursor_base
%6.8.19

%Debug check:
clc
b2base_num_samples  = 2; 
cursor              = [0 1 1 0 1 0 0.1 0 1 0 0 1]; 
base_val            = 0.5; 
hit_idxs            = find(cursor >= base_val)

[valid_hit_idxs] = ...
    analyze_hits_from_cursor_base(cursor, hit_idxs, base_val, b2base_num_samples);
valid_hit_idxs