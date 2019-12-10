%bin_ND_tb

%1) pick the bin edges

x_bin_edges = [0.5 1.5 2.5];
y_bin_edges = [3.5 4.5 5.5 6.5];

num_x_bins = length(x_bin_edges)+1;
num_y_bins = length(y_bin_edges)+1;

%2) pick a test histogram, generate data for the histogram.

truth_hist = zeros(num_y_bins, num_x_bins);