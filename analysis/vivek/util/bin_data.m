function [data_binned, bin_counts] = bin_data(bin_edges, data)
%input:
%data = vector of size n x 1
%bin_edges = vector of size 1 x num_bins
%output:
% data_binned = vector of size n x 1
% bin_counts = vector of size = num_bins x 1
% 
% num_bins = length(bin_edges)+1 
% minimum possible value of data_binned is 1
% maximum possible value of data_binned is num_edges + 1;
%
num_bins = (length(bin_edges)-1)+2;
data_binned = zeros(size(data));
bin_counts = zeros(num_bins,1);

for i=1:length(data)
    data_binned(i) = bin_value(bin_edges, data(i))+1;
end

for i=1:num_bins
    bin_counts(i) = sum(data_binned == i);
end