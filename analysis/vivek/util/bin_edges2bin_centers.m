function [bin_centers] = bin_edges2bin_centers(bin_edges)
%input:
%bin_edges = vector 1 x n
mean_bin_width = mean(conv(bin_edges, [1 -1], 'valid'));
bin_centers = [bin_edges(1)-mean_bin_width/2 .5*conv(bin_edges, [1 1], 'valid')  bin_edges(end)+mean_bin_width/2];