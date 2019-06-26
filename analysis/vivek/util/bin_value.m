function result = bin_value(bin_edges, value)
%Convention:
%num_bins = length(bin_edges)+1 
%if value is less than the lowest bin_edge, then its bin is 0
%if value is greater than largest bin_edge, then its bin is
%length(bin_edges)
 

cmp = value > bin_edges;
if(isempty(cmp))
    result = 0;
else
    result = sum(cmp);
end