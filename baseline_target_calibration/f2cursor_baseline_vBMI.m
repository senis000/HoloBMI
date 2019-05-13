% function [cursor_obs] = f2cursor_baseline_vBMI(f_data, bData, prefix_win, f0_win, dff_win)
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
% baselineDataFile = 'BaselineOnline190512T025909.mat'
baselineDataFile = 'BaselineOnline190512T022047.mat'

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


% % FROM BMI:
% if data.frame == initFrameBase
%     disp('blerg')
%     baseval = single(ones(numberNeurons,1)).*unitVals;
% elseif data.frame <= baseFrames
%     baseval = (baseval*(data.frame - 1) + signal)./data.frame;
% else
%     baseval = (baseval*(baseFrames - 1) + signal)./baseFrames;
% end


%%
%Calculate f0 as in BMI: 
num_samples = size(f_raw,1)
f0 = zeros(num_samples,num_neurons); 
for i=1:length(f0)
    if i==1
        f0(i,:) = f_raw(i,:);
    elseif i<f0_win
        f0(i,:) = (f0(i-1,:)*(i-1)+f_raw(i,:))/i;
    else
        f0(i,:) = (f0(i-1,:)*(f0_win-1)+f_raw(i,:))/f0_win;
    end
end
%Truncate data based on the f0_win:
f_postf0 = f_raw;
f0_mean = repmat(nanmean(f_postf0, 1), size(f_postf0,1), 1);

%%
%Second, smooth::
num_samples = length(f_postf0); %max(length(dff_z)-(dff_win-1), 0); 
f_smooth = zeros(num_samples, num_neurons); 
smooth_filt = ones(dff_win,1)/dff_win; 
for i=1:num_neurons
    f_smooth(:,i) = conv(f_postf0(:,i), smooth_filt, 'same'); 
end

%%
%Second, compute dff and dff_z:
dff_smooth = (f_smooth-f0)./f0;
%mean center the dff:
%n_mean = nanmean(dff,1); %1 x num_neurons
mean_mat = repmat(bData.n_mean, size(dff_smooth,1), 1);
dffc = dff_smooth-mean_mat;
%divide by std:
% n_std = nanstd(dffc, 0, 1); %var(dffc, 0, 1).^(1/2); %1 x num_neurons
dff_z = dffc./repmat(bData.n_std, [size(dff_smooth,1) 1]);

n_analyze = dff_z; 

%%
first_samples = 1:f0_win;
last_samples = (f0_win+1):length(f_smooth); 

n_first = dff_smooth(first_samples,:); 
n_last = dff_smooth(last_samples,:); 

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


%%
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
% n_mean =
% 
%   Columns 1 through 4
% 
%    -0.0251    0.0045    0.0450    0.0365
% 
%   Columns 5 through 8
% 
%     0.0130    0.0344    0.0282    0.0164
% 
% 
% n_std =
% 
%   Columns 1 through 4
% 
%     0.4837    0.6212    0.7363    0.9009
% 
%   Columns 5 through 8
% 
%     0.3374    0.5126    0.1443    0.4741



%%
%Cursor: 
cursor_obs = n_analyze*bData.decoder; 

h = figure;
plot(cursor_obs); 
title('Cursor from base BMI'); 