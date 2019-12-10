function [data_binned, grid_bin_counts] = bin_ND(data, edges_cell)
%4.29.14
%input:
%data           - N x D, N = number of data points, D = dimensionality of data
%edges_cell     - 1 x D, cell array.  contains bin edges for the particular
%dimension.
%
%output:
%data_binned    - N x D, each entry has a bin value
%grid_bin_counts     - size: (num_edges_D1+1) x (num_edges_D2+1)
%ex:
%grid_bin_counts(2,2) = num_counts between edges_cell{1}(1:2) and
%edges_cell{2}(1:2)
%
% binning:
% i specify num_edges.  num_bins = num_edges + 1.  the base case is that
% you specified one edge, which yields two bins on either side of the edge.



N = size(data,1);
D = size(data,2);

num_bins_vec    = zeros(1,D);
for i=1:D
    %NUM_BINS = NUM_EDGES + 1, by this code's convention.
    num_bins_vec(i)    = length(edges_cell{i})+1;
end

grid_bin_counts     = zeros(num_bins_vec);

data_binned = zeros(N,D);
%function [data_binned, bin_counts] = bin_data(bin_edges, data)
for i=1:D
    edges_i                     = edges_cell{i}; %vector of edges for the binning the ith dimension
    [data_binned_i, bin_counts] = bin_data(edges_i, data(:,i));
    data_binned(:,i)            = data_binned_i;
end

% N;
% D;
% size(data_binned);
ind = submat2ind(num_bins_vec, data_binned);

for i=1:N
    grid_bin_counts(ind(i)) = grid_bin_counts(ind(i)) + 1;
end