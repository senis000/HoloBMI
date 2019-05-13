
%%
%BASE:
% n_mean =
% 
%   Columns 1 through 4
% 
%    -0.0055   -0.3439   -0.2651   -0.3972
% 
%   Columns 5 through 8
% 
%    -0.6538   -0.5267   -0.6554   -0.5452
% 
% 
% n_std =
% 
%   Columns 1 through 4
% 
%     0.4908    0.4092    0.5184    0.5168
% 
%   Columns 5 through 8
% 
%     0.1252    0.2871    0.0565    0.2200

%%
%BASE BMI:
n_mean = nanmean(dff_smooth,1) %1 x num_neurons
n_std = nanstd(dff_smooth, 0, 1)
%zscore values on the smoothed dff:
% n_mean =
% 
%   Columns 1 through 4
% 
%    -0.0251    0.0046    0.0451    0.0367
% 
%   Columns 5 through 8
% 
%     0.0130    0.0345    0.0282    0.0164
% 
% 
% n_std =
% 
%   Columns 1 through 4
% 
%     0.4792    0.6176    0.7322    0.8964
% 
%   Columns 5 through 8
% 
%     0.3308    0.5053    0.1272    0.4692

%%
n_mean = nanmean(n_first,1) %1 x num_neurons
n_std = nanstd(n_first, 0, 1)

% n_mean =
% 
%   Columns 1 through 4
% 
%    -0.1506    0.0130    0.0583    0.2665
% 
%   Columns 5 through 8
% 
%    -0.0104    0.0166    0.0019    0.0267
% 
% 
% n_std =
% 
%   Columns 1 through 4
% 
%     0.4108    0.6429    0.6563    1.0746
% 
%   Columns 5 through 8
% 
%     0.2729    0.1792    0.0916    0.4615

%%
n_mean = nanmean(n_last,1) %1 x num_neurons
n_std = nanstd(n_last, 0, 1)

% n_mean =
% 
%   Columns 1 through 4
% 
%    -0.0054    0.0033    0.0431    0.0007
% 
%   Columns 5 through 8
% 
%     0.0167    0.0373    0.0323    0.0148
% 
% 
% n_std =
% 
%   Columns 1 through 4
% 
%     0.4861    0.6135    0.7434    0.8597
% 
%   Columns 5 through 8
% 
%     0.3388    0.5387    0.1314    0.4704



