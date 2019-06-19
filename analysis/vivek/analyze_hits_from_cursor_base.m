function [valid_hit_idxs] = ...
    analyze_hits_from_cursor_base(cursor, hit_idxs, base_val, b2base_num_samples)
%6.8.19
%posthoc analyzes hit times, and only returns the hit times that are after
%a back2base event
%b2base_event logic: 
%if between hits the number of samples below base_val is >=
%b2base_num_samples

valid_hit_idxs = [hit_idxs(1)]; 
for i=2:length(hit_idxs)
    hit_interval    = hit_idxs(i)-hit_idxs(i-1)-1;
    if(hit_interval >= b2base_num_samples)
        interval_idxs   = hit_idxs(i-1):hit_idxs(i);
        cursor_interval = cursor(interval_idxs);         
        
        b2base_bool = sum(cursor_interval <= base_val) >= b2base_num_samples; 
        if(b2base_bool)
            valid_hit_idxs = [valid_hit_idxs hit_idxs(i)]; 
        end
    end
end