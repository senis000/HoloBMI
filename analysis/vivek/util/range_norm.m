function [data_range_norm, data_min, data_max, data_range] = range_norm(data, max_percentile, min_percentile)
%To Do: extend to multidimensional 
data_min = prctile(data, min_percentile); 
data_max = prctile(data, max_percentile); 
data_range = data_max - data_min; 
data_range_norm = (data-data_min)/data_range; 