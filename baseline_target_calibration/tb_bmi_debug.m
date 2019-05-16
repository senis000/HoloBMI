
%5.14.19
%%
% function BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, cursor_zscore_bool, debug_bool)

%%
debug_bool = 1;

frameRate = 30
expt_str = 'BMI'
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

vectorHolo = []
vectorVTA = []

folder = 'E:\VivekNuria\expt\HoloBmi\debug'
% folder = '/Users/vivekathalye/Dropbox/Data/holo_bmi_debug/190513'
animal = 'N3'
day = 'test'

cursor_zscore_bool = 0; 

%BMI Calibration File:
baselineCalibrationFile = 'BMI_target_info_20190515T152926.mat'
savePath = fullfile(folder, animal, day); %[folder, animal, '/',  day, '/'];
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
exist(fullfile(savePath, baselineCalibrationFile))
%%
test = load(fullfile(savePath, baselineCalibrationFile))

%%
%Debug Input File with fluorescence: 
bmi_file = fullfile(savePath, 'BMI_online190514T021302.mat'); 
exist(bmi_file)

bmi_data = load(bmi_file); 
F = bmi_data.data.bmiAct; %num_neurons X num_samples
%1) remove nans
F(:,isnan(F(1,:))) = []; 
%2) pad front with baseFrames
baseFrames = 40; 
F_pad = zeros(size(F,1), baseFrames); 
F = [F_pad F]; 

debug_input_file = fullfile(savePath, 'debug_input.mat'); 
save(debug_input_file, 'F'); 
test = load(debug_input_file)
%%
num_E = 6;
num_E1 = 3; 
num_E2 = 3; 

n_prefix = 100*ones(num_E,40); 
n_base = 100*ones(num_E,100); 
n_T = [zeros(num_E1,5); 10000*ones(num_E2,5)];
n_z = zeros(num_E,10); 

F_v = [n_prefix n_base n_z n_T n_z n_T n_z n_T n_z n_z n_z];

%%
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
offset = 0; 
[h, offset_vec] = plot_E_activity(F_v', [ones(4,1); 2*ones(4,1)], E_color, offset)

%%
% baseSeed = ones(num_E, 1)+nan; 
baseSeed = 10000*ones(num_E, 1); 

%% 
BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, expt_str, ...
    baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, cursor_zscore_bool, debug_bool, F_v, baseSeed)


% function BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, 
%expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, cursor_zscore_bool, debug_bool, debug_input)

%%
%Plot base: 
base_valid = data.baseVector(:, ~isnan(data.baseVector(1,:))); 
h = figure;
plot(base_valid(1,:))

%%
debug_bool = 1;

frameRate = 30
expt_str = 'HoloVTA_pretrain'
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

vectorHolo = []
vectorVTA = []

% folder = 'F:\VivekNuria\expt\HoloBmi'
folder = '/Users/vivekathalye/Dropbox/Data/holo_bmi_debug/190513'
animal = 'NY20'
day = 'test'

cursor_zscore_bool = 0; 

%BMI Calibration File:
baselineCalibrationFile = 'BMI_target_info_20190514T015556.mat'
savePath = fullfile(folder, animal, day); %[folder, animal, '/',  day, '/'];
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
exist(fullfile(savePath, baselineCalibrationFile))

%%
%Debug Input File with fluorescence: 
bmi_file = fullfile(savePath, 'BMI_online190514T021302.mat'); 
exist(bmi_file)

bmi_data = load(bmi_file); 
F = bmi_data.data.bmiAct; %num_neurons X num_samples
%1) remove nans
F(:,isnan(F(1,:))) = []; 
%2) pad front with baseFrames
baseFrames = 40; 
F_pad = zeros(size(F,1), baseFrames); 
F = [F_pad F]; 

debug_input_file = fullfile(savePath, 'debug_input.mat'); 
save(debug_input_file, 'F'); 
test = load(debug_input_file)
%%
n_prefix = 100*ones(8,1); 
n_prefix_len = repmat(n_prefix, 1, 40); 

n_base = 100*ones(8,1); 
n_base_len = repmat(n_base, 1, 100);

n_T = [zeros(4,1); 10000*ones(4,1)];
n_z = zeros(8,1); 
n_E1 = [10000*ones(4,1); zeros(4,1)]; 
n_E1_len = repmat(n_E1, 1, 10); 

F_v = [n_prefix_len n_base_len repmat(n_z, 1, 10) n_E1_len repmat(n_z, 1, 1) n_T repmat(n_z, 1, 30) repmat(n_T, 1, 20)]; 

vectorHolo = [length(n_prefix_len) + ...
    length(n_base_len) + ...
    10 + ...
    1 170]; 

%%
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
offset = 0; 
[h, offset_vec] = plot_E_activity(F_v', [ones(4,1); 2*ones(4,1)], E_color, offset)
vline(vectorHolo); 

%% 
expt_str = 'BMI'
BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, expt_str, ...
    baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, cursor_zscore_bool, debug_bool, F_v)

%%
result = load(fullfile(savePath, 'BMI_online190514T200139.mat')); 

%%
[h, offset_vec] = plot_E_activity(F_v', [ones(4,1); 2*ones(4,1)], E_color, offset)
vline(result.data.vectorHoloCL); 


