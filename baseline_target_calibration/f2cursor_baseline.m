function cursor_obs = f2cursor_baseline(f_data, bData, prefix_win, f0_win, dff_win)
%5.12.19
%Function which converts neural fluorescence to a cursor value, using the
% algorithm from the baseline calibration
%

%f_data: num_neurons X num_time_samples
%Remove nan:

%%
% Debug Inputs:
% folder = 'F:\VivekNuria\expt\HoloBmi'
frameRate = 30
expt_str = 'BMI'
movingAverageFrames = 2

folder = '/Users/vivekathalye/Dropbox/Data/holo_bmi_debug_190512'
animal = 'test'
day = 'test'

baselineCalibrationFile = 'BMI_target_info.mat'
baselineDataFile = 'BaselineOnline190512T022047.mat'
% baselineDataFile = 'BaselineOnline190512T025909.mat'

savePath = fullfile(folder, animal, day); %[folder, animal, '/',  day, '/'];
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
exist(fullfile(savePath, baselineCalibrationFile))
exist(fullfile(savePath, baselineDataFile))

bData = load(fullfile(savePath, baselineCalibrationFile)); 
load(fullfile(savePath, baselineDataFile)); 
%baseActivity
%removenan:
f_data = baseActivity; 
f_data(:,isnan(f_data(1,:))) = [];
num_samples = size(f_data, 2); 


prefix_win = 40;
f0_win = round(2*60 * frameRate);
dff_win = movingAverageFrames;

%%
f_data(:,isnan(f_data(1,:))) = [];
f_base = f_data.';

%Throw out prefix frames:
E1_raw = f_base((prefix_win+1):end,bData.E1_base); 
E2_raw = f_base((prefix_win+1):end,bData.E2_base); 
f_raw = [E1_raw E2_raw]; %first E1, then E2

%2) Ensemble information
num_E1 = length(bData.E1_base); 
num_E2 = length(bData.E2_base); 
num_neurons = num_E1 + num_E2;

E_id = [1*ones(num_E1, 1); 2*ones(num_E2, 1)]; 
E1_sel = E_id==1; 
E1_sel_idxs = find(E1_sel); 
E2_sel = E_id==2; 
E2_sel_idxs = find(E2_sel); 


%%
%Calculate f0:

%Implemented in the same fashion as in the BMI currently
num_samples = size(f_raw,1);
f0 = zeros(num_samples-f0_win+1, num_neurons); 
f0(1,:) = mean(f_raw(1:f0_win, :), 1);
for i = 2:length(f0)
    f0(i,:) = f0(i-1)*((f0_win-1)/f0_win) + f_raw((i+f0_win-1), :)/f0_win; 
end
%Truncate data based on the f0_win:
f_postf0 = f_raw(f0_win:end, :); 
f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);

%%
%Second, compute dff and dff_z:
dff = (f_postf0-f0)./f0;
%mean center the dff:
% n_mean = nanmean(dff,1); %1 x num_neurons
mean_mat = repmat(bData.n_mean, size(dff,1), 1);
dffc = dff-mean_mat;
%divide by std:
% n_std = nanstd(dffc, 0, 1); %var(dffc, 0, 1).^(1/2); %1 x num_neurons
dff_z = dffc./repmat(bData.n_std, [size(dff,1) 1]);

%%
n_mean = nanmean(dff,1) %1 x num_neurons
% bData.n_mean

n_std = nanstd(dff,0,1) %1 x num_neurons
% bData.n_std

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
%Third (final), smooth dff:
num_samples = max(length(dff_z)-(dff_win-1), 0); 
n_analyze = zeros(num_samples, num_neurons); 
smooth_filt = ones(dff_win,1)/dff_win; 
for i=1:num_neurons
    n_analyze(:,i) = conv(dff_z(:,i), smooth_filt, 'valid'); 
end

%%
%Cursor: 
cursor_obs = n_analyze*bData.decoder; 

%%
h = figure;
plot(cursor_obs); 

