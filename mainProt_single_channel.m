%%%TODO'S!!!!
%{
Base: disp number of target hits (without b2base)
BMI: disp number of delivered holo stims
streamline analysis of baseline, pretrain, bmi data, plots in same format
rewards per minute plot
anticipatory licking?
Baseline CA vs stim time
 Close-loop E2 not being high-ish
analyze baseline E2
 what window around stim triggers reward
%
%For Experiment:
%A protocol for determing reliably stimmed cells, and determing power and
duration of stim
%}

%% Main protocol for the experiment
%--------------------------------------------------------------------------
%BEFORE ANIMAL IN BOX:
%DO:
% Hook up BNCs: 
% 1) BMI solenoid, AI5
% 2) Monaco Trig, AI6
% 3) Frame Trig, AI7
% 4) Holo Trig PFI1
% 
% Power Arduino: 
%   (Power supply needed to power solenoid, can't control solenoid on USB power)
% Voltage Recording: All Inputs Active (check 6+7)
%
%Fill syringe with sucrose cuz of gravity
%Run test_bmi_nidaq_triggers.m
%   check triggers for 1) get image, 2) trig photostim, 3) trig reward
%calibrate solenoid opening time
%
%load pyctrl expt for the mouse
%collect load cell sensor baseline data
%
%Put mouse in
%put gel from headbar to ear
%adjust spout so mouse can lick
%put objective
%--------------------------------------------------------------------------
home_dir = 'G:\VivekNuria\Code\HoloBMI'
cd(home_dir)
addpath('G:\VivekNuria\Code\HoloBMI\roi_process'); 
addpath('G:\VivekNuria\Code\HoloBMI\baseline_target_calibration'); 
addpath('G:\VivekNuria\Code\HoloBMI\matlab_exchange'); 
addpath(genpath('G:\VivekNuria\Code\dex')); 
addpath(genpath('G:\VivekNuria\Code\Kenichi')); 


% define Animal, day and folder where to save
animal = 'NY35'; day = 'D1_test';
folder = 'E:\vivek\190907';

% define posz TODO can we get this from prairie?
posz = 0;

pl = actxserver('PrairieLink.Application');
pl.Connect(); disp('Connecting to prairie');
% Prairie variables
px = pl.PixelsPerLine();
py = pl.LinesPerFrame();

micronsPerPixel.x = str2double(pl.GetState('micronsPerPixel', 'XAxis')); 
micronsPerPixel.y = str2double(pl.GetState('micronsPerPixel', 'YAxis')); 

pl.Disconnect();

% px = 512; 
% py = 512;
frameRate = 29.989;

chan_data = struct(...
    'label', 'g', ...
    'chan_idx', 2); %in RGB, G is 2nd

onacid_bool = false;

savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%%
%--------------------------------------------------------------------------
%DO:
%find FOV
%-disable the motor control!!!!
%-set imaging FOV.  Use High Power to calc image for holoMask
%-lower imaging power back to normal
%TODO: have the script take a 1000 frame video and calculate the mean
%--------------------------------------------------------------------------
%% Select gcamp neurons
%--------------------------------------------------------------------------
%DO:
%-mainMasksPrairie, make holoMask
%-deleteMask, delete bad neurons.  
%-addcell, add neurons
%-add green cells: 
%--[holoMaskRedGreen, ~,~] = addcell (Img, holoMask,9);
%--[holoMaskRedGreen, ~,~] = addcell (Img, holoMaskRedGreen,9);
%--------------------------------------------------------------------------
im_summary_path    = fullfile('G:\vivek\190822_NY35_good_stim_tests\NY35\D1_test', 'green_std.tif'); 
% fullfile('E:\vivek\190822\NY35\D1_test', 'chan_mean.tif'); 
exist(im_summary_path)
im_summary = imread(im_summary_path); 
% holoMask = zeros(size(green_im)); 

%%
%Scale the image, in order to help see ROIs better.
%If no modification to original image needed, just run code, in 
% 'scale_im_interactive' set min_perc = 0, max_perc = 100
im_sc_struct = struct(...
    'im', [], ...
    'minmax_perc', [], ...
    'minmax', [], ...
    'min', [], ...
    'min_perc', [], ...
    'max', [], ...
    'max_perc', []); 
num_im_sc = 0; 
[im_sc_struct, num_im_sc] = scale_im_interactive(im_summary, im_sc_struct, num_im_sc);

%%
%INITIALIZE ROI DATA
%--------------------------------------------------------------------------
%DO:
%find FOV
%-choose im_bg, set the index of 'im_sc_struct' 
%--------------------------------------------------------------------------

%CHOOSE IM_BG
im_bg = im_sc_struct(end).im; 

h = figure;
imagesc(im_bg); 
axis square
colormap('gray'); 
title('selected background image for identifying ROI'); 

%plot_images contains a set of images so user can tell if ROI selection is
%appropriate.
plot_images = struct('im', [], 'label', ''); 
plot_images(1).im = im_summary; 
plot_images(1).label = 'green std';

plot_images(2).im = im_bg; 
plot_images(2).label = 'scaled';

%%
%INIT ROI_DATA
auto_init = 1;  %initializes roi_data using automatic cell detection: 
%Parameters for auto cell detection:
%Following were for zoom=2 on bruker:
template_diam = 25; %diamter of difference of Gaussians in pixels
thres = 0.5; %cell detection threshold as correlation coefficient
cell_diam = 7; %CELL_DIAM is diameter used for dilation.
finemode = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
temmode = 0; 
if auto_init
    %FIND ROI AUTOMATICALLY 
    [mask_intermediate, ~] = imFindCellsTM (im_bg, template_diam, thres, cell_diam, finemode, temmode);
    init_roi_mask = bwlabel(mask_intermediate);
    findCenter (init_roi_mask, im_bg);
    roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, chan_data);
else
    roi_data = init_roi_data(im_bg); 
end

%%
%FOR TESTING: 
roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, chan_data);

%%
%FOR TESTING:
[roi_ctr] = roi_bin_cell2center_radius(roi_data.roi_bin_cell);

%%


%%
%FOR TESTING: 
screen_size = get(0,'ScreenSize');
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
% hold on;
imagesc(roi_data.im_roi); %colormap('gray');  
axis square;
% scatter(roi_data.x, roi_data.y, pi*roi_data.r.^2, 'r'); 

%%
%Delete ROI if needed: 
disp('Deleting ROIs from image!');
[roi_data] = delete_roi_2chan(plot_images, roi_data);

%%
disp('Adding ROIs to image!'); 
[roi_data] = draw_roi_g_chan(plot_images, roi_data);

%%
%See the holoMask:
holoMask = roi_data.roi_mask; 
h = figure;
imshow(holoMask)

%% Save roi_data
%TODO: 

%%
filetosave = fullfile(savePath, 'red.mat');
load(filetosave)
h = figure;
imshow(holoMask)

%% prepare HOLO STIM of individual neurons
%Create gpl, xml files for individual points

% load environment
env_path = fullfile('G:\VivekNuria\utils', "Tseries_VivekNuria_holo_all.env")
exist(env_path)


%%
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%%
%Darcy recommends 5-10 spirals
% size_um = roi_data.
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', 100, ...
    'SpiralSize', 0.7, ...
    'SpiralRevolutions', 10); 
markpoints_data = repmat(init_markpoints, [roi_data.num_rois 1]); 
for roi_i = 1:roi_data.num_rois
    spiral_size_um  = 2*roi_data.r(roi_i)*micronsPerPixel.x;
    spiral_size     = spiral_size_conversion*spiral_size_um;
    markpoints_data(roi_i).SpiralSize = spiral_size; %double(ceil(10*spiral_size)/10);
end
%%
% creates holos
createGplFile_v2(savePath, markpoints_data, roi_data.x, roi_data.y, posz, roi_data.r, px)

%define where to save the file
%TODO: fix this code
savePrairieFilesHolo(savePath)

% create the stim train for prairie
num_neurons_stim = roi_data.num_rois; %max(max(holoMask));

numberNeurons=max(max(holoMask));
power = 0.2; 
duration = 30; 
numSpirals = 10; 

initDelay = 2000; %(ms) time bw stim delivery
iterations = 1;
reps = 1;

powerVector = power*ones(1,numberNeurons); %*0.1; %0.2 = 50, 0.1=25
durationVector = duration*ones(1,numberNeurons);
spiralVector = numSpirals*ones(1,numberNeurons);

createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector, iterations, '', false)
% createXmlFile(savePath, numberNeurons, reps, power, varName, flagRandom)

% Update prairie view repetitions based on num neurons to stim
num_stim_neurons = max(max(holoMask));
stim_time_per_neuron = 2.1;
disp(num_stim_neurons*stim_time_per_neuron*frameRate)

%--------------------------------------------------------------------------
%DO: 
% upload .gpl in MarkPoints (Top half)
% upload .xml in MarkPoints (Bot half)
% update T-series repetitions in Prairie View with above number
% Make sure Voltage Recording has all channels enabled
%--------------------------------------------------------------------------

%% Run HOLO STIM to check stim-able neurons
%This stims one neuron at a time.
%--------------------------------------------------------------------------
%DO: 
%1) Do live scan check
%--------------------------------------------------------------------------

clear s
HoloAcqnvsPrairie(folder, animal, day, holoMask)
%WE WANT RED+GREEN FOR THIS
%TODO: make this closed loop, and wait for the neurons to be inactive
%before stimming them.

%--------------------------------------------------------------------------
%DO: 
%Image-block ripping utility: Convert the holostim acqn
%Load holoOnlineX.mat and voltage recording for plotting
%(Import csv to matlab: output type is Numeric Matrix.  Name it "voltageRec")
%--------------------------------------------------------------------------

%% Obtain spatial components
% run OnAcidPrairieNotebook.ipynb 
%TODO: We still have to confirm if onacid on holostim gives same spatial 
%components as onacid on baseline

%COPY THE FOLLOWING LINES INTO ANACONDA:

% obtain A from onacid, compare to red neurons and bring it to matlab
% in python run OnAcid_Prairie_holo and obtain_components
% "obtain_componenents.py" -> 
%         'AComp' : sparse matrix with the spatial filters,
%         'CComp' : temporal filters
%         'com' : position of neurons,
%         'redlabel': labels with true on red neurons,
%         'redIm' : image of the red channel,
%         'baseIm' : background of the image given by caiman
% it also saves figures for sanity check

% while onacid does its magic 

%%
%(Image-Block Ripping Utility) Convert holostim file with bruker converter 
% load the VoltageRec to check the results of holoStim
plotHoloStimTimeLock(holoActivity, voltageRec, 40, 1000) % --> To plot the result of
%ToDo: for plotting, do sliding window deltaf/f
%ToDo: allow us to select the idxs of neurons to plot
%ToDo: a version that just plots each individual neuron, we type 'Y' or 'N'
%to make it a candidate

%%
%--------------------------------------------------------------------------
%DO: 
%Select E2_base, choose neurons which appear to be stimmed well.  You can
%choose more than 4 neurons, and then re-run once you've chosen your 4.  
%--------------------------------------------------------------------------
% E2_base = sort([21    36   127   196], 'ascend')
% E2_base = sort([58 67 56 27 87 89 24 51 8 95 72 91 16 115 98 75 153], 'ascend')

% E2_base = ...
%     [8 2 7 14 20 23 16];

% E2_base = E2_base([3 5 7])

E2_base = [7 20 16]; 
 
%  E2_base = [7     6    47    24]; 
% E2_base = sort([7     3    47    57], 'ascend')
% E2_base = [3 7 47 57]; 
% E2_base = [98   158   161   164]; 
% E2_base = sort([14 10 26 9], 'ascend')

%TODO: update protocol to iterate from E2_candidates to E2_base

%% Holo stim of Ensemble (only stim cells)
%TODO create file + creating the xml etc
% create the masks for holo stim the 4 ensemble neurons
if onacid_bool
    StimMask=deleteNonEnsemble (AComp, E2_base, px, py);
else
    StimMask = holoMaskRedGreen;
    StimMask(~ismember(StimMask,E2_base))= 0;
end
figure;
imshow(StimMask); 

%%
%GPL for Stim Ensemble
pl = actxserver('PrairieLink.Application');
pl.Connect();
savePrairieFiles(savePath, pl, '4neurons_together')
createGplFile(savePath, StimMask, posz, px, 'ensemble_')

createBot(savePath, x(E2_base),y(E2_base))

loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);

%laser power 50, 4 iterations with
% 10 sec bw
% CHANGE NUMBER OF REPS IN T-SERIES:
%4 iterations of 10 sec = 4*10*30
num_iter = 4
time_bw_iter = 10
imaging_reps = num_iter*time_bw_iter*30


%--------------------------------------------------------------------------
%DO: 
%1) upload the GPL file
%2) add points to Point series
%3) Click BOT in Image Window, load BOT.cfg, only display ROI of interest
%4) Run BOT, and adjust (duration, power) of each neuron's stim
%--------------------------------------------------------------------------
%TODO: a better automated way of determining the stim-able neurons

%% Power/Duration for each candidate neuron
%Scratch to save the power/duration parameters.
%{
like:

1
NO

2
NO

3
GOOD
p=15
d=20

4
NO

5
yes
p=15
d=20

6
NO

7
YES
p=35
d=20

%}
%%
% E2_base = E2_base([2 11 9 15])

% 3
% GOOD
% p=15
% d=20
% 
% 5
% yes
% p=15
% d=20
% 
% 7
% YES
% p=35
% d=20

%%
powerVector = .004*[15 15 35];
durationVector = [20 20 20];

%TODO: save the data of the single cell stim!!
%%
s = daq.createSession('ni');
addDigitalChannel(s,'dev5','Port0/Line0:2','OutputOnly');
ni_out = [0 0 0]; 
outputSingleScan(s,ni_out);%set   
ni_getimage = [1 0 0]; 
ni_reward   = [0 1 0]; 
ni_holo     = [0 0 1]; 
%%
outputSingleScan(s,ni_holo);
pause(0.01); outputSingleScan(s,ni_out)

%% Run 4 masks together (obsolete)
% pl.SendScriptCommands("-ts"); 
% pause(1);

%%
clear s
%% Baseline acquisition
%Note: loads the result of OnAcid / holoMask
%Do this after we confirm we can stim some cells

%--------------------------------------------------------------------------
%DO: 
%Remove Red Channel from Image Window 1 (prairie view).
%0) (zero pmt+power) put water
% check FOV didn't move
%1) start video
%2) start load cells
%3) start pyctrl
%4) Run this cell
%--------------------------------------------------------------------------

if ~onacid_bool
    AComp = 0;
else
    load(fullfile(savePath,'redcomp.mat'));
end


% Baseline environment already removes MARKPOINTS and set the reps to 27000
BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMask, onacid_bool, frameRate);
% BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMaskRedGreen, onacid_bool, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

%--------------------------------------------------------------------------
%D0:
%0) Abort T-series (cuz of voltage recording)
%1) pyctrl stop
%2) load cells stop
%3) video stop
%4) Drag load cell data to folder
%5) Drag video to folder
%--------------------------------------------------------------------------

%% Selection of neurons
% plots neurons so we can select which ones we like the most 

%--------------------------------------------------------------------------
%D0:
%1) Copy paste base_file name 'BaselineOnlineX.mat'
%--------------------------------------------------------------------------

% load by hand! --> (you can blame Vivek for this :P load(fullfile(savePath,'BaselineOnline.mat'));
if onacid_bool
    totalneurons = min(size(AComp,2), 20);
else
    totalneurons = max(max(holoMaskRedGreen));
    CComp = [];
    YrA = []; 
end

% totalneurons = 40; 
%Copy paste base_file path: 
base_file = fullfile(savePath, 'BaselineOnline190526T113422.mat')
load(base_file); 
% plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)
plotNeuronsBaseline(baseActivity, CComp, YrA, 30)

%TODO:  
%ToDo: for plotting, do sliding window deltaf/f
%ToDo: allow us to select the idxs of neurons to plot (so we can show the
%neurons we chose for the BMI)
%ToDo: a version that just plots each individual neuron, we type 'Y' or 'N'
%to make it a candidate
%--------------------------------------------------------------------------
%D0:
%1) Choose E1_base
%--------------------------------------------------------------------------
%%
%Manually enter:
E1_base = sort([21 11 18], 'ascend') %JUST NEEDS GCAMP
% E2_base = sort([35 37 13 4], 'ascend')

% E1_base = sort([64 19 141 83], 'ascend') %JUST NEEDS GCAMP
% E2_base = sort([7 12 34 48], 'ascend')

ensembleNeurons = [E1_base, E2_base];
plotNeuronsEnsemble(baseActivity, ensembleNeurons)

% E2_candidates = [39 45 59 37 88 6 26 46 78 48 22 20 33]

%%
%OPTION: Use previously collected BMI data as the baseline data: 
%
% bmi_file = fullfile(savePath, 'BMI_online190515T010526.mat'); 
% bmi_data = load(bmi_file); 
% bmi_base = fullfile(savePath, ['base_' 'BMI_online190515T010526.mat']);
% baseActivity = bmi_data.data.bmiAct(:, ~isnan(bmi_data.data.bmiAct(1,:))); 
% save(bmi_base, 'baseActivity'); 
% 
% E1_base = [1 2 3 4]; 
% E2_base = [5 6 7 8]; 

%% Calibrate Target with Baseline simulation
%--------------------------------------------------------------------------
%D0:
%1) Parameters: 
% - sec_per_reward_range
% - f0_win (F0: how many frames to average over)
% - dff_win (F for Dff: how many frames to average over)
%--------------------------------------------------------------------------

% base_file = fullfile(savePath, 'BaselineOnline190514T221822.mat')
% base_file = bmi_base; 

exist(base_file)
n_f_file = base_file;
exist(n_f_file)
ndata = load(n_f_file);
num_base_samples = sum(~isnan(ndata.baseActivity(1,:))); 
baseline_frameRate = num_base_samples/(15*60);
A_file = fullfile(savePath, 'red.mat'); 
exist(A_file)
onacid_bool = 0

sec_per_reward_range = [120 80]; 

frames_per_reward_range = sec_per_reward_range*baseline_frameRate %[1 1.5]*60*frameRate
%multiply by frames per minute to convert
%to 

target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 2*60*ceil(frameRate)
dff_win_bool = 1
dff_win = 4
 
reward_per_frame_range = 1./frames_per_reward_range

cursor_zscore_bool = 0;
f0_init_slide = 0; 

close all
%  baseline2target_vBMI(n_f_file, A_file, onacid_bool,  ...
%     [1 2 3 4], [5 6 7 8], frames_per_reward_range, target_on_cov_bool, ...
%     prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
%     cursor_zscore_bool, f0_init_slide);



%  baseline2target_vBMI(n_f_file, A_file, onacid_bool,  ...
%     E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
%     prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
%     cursor_zscore_bool, f0_init_slide);

baseline2target_vE1strict(n_f_file, A_file, onacid_bool,  ...
    E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
    prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
    cursor_zscore_bool, f0_init_slide)

%ToDo: return the filename

% baseline2target(n_f_file, A_file, onacid_bool, E1_base, E2_base, frames_per_reward_range, ...
%     target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath)


% run the simulation of baseline
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 
%To Do: Show the percent correct of the pretrain period, based on the
%calibration. 

%--------------------------------------------------------------------------
%D0:
%Note down: 
% - T value
% T: 0.54
% num_valid_hits: 7
% num_hits: 74
%--------------------------------------------------------------------------

% %%
% t = load(fullfile(savePath, 'target_calibration_ALL_20190524T122903.mat'))

if onacid_bool
    StimMask=deleteNonEnsemble (AComp, E2_base, px, py);
else
    StimMask = holoMaskRedGreen;
    StimMask(~ismember(StimMask,E2_base))= 0;
end

figure;
imshow(StimMask); 

pl = actxserver('PrairieLink.Application');
pl.Connect();
savePrairieFiles(savePath, pl, '4neurons_together')
createGplFile(savePath, StimMask, posz, px, 'ensemble_')

createBot(savePath, x(E2_base),y(E2_base))

loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);

%ACTIONS:
% SELECT BOT on image window
% upload the GPL file
% upload the BOT file
% Manually sselect neurons, make group, laser power 50, 4 iterations with
% 10 sec bw
% CHANGE NUMBER OF REPS IN T-SERIES
%% Run 4 masks together

pl.SendScriptCommands("-ts"); 
pause(1);

%%
%Power/Duration Reminder: 
%
%idx1: i_base = 7
% p15,d20
%
% id2: i_base= 12
%p10, d10
%
% id3, i_base = 34
%p10, d20

% id3, i_base = 48
%p10, d10
%ToDo: automatically convert to power setting
%% Holo stim checking connectivity
% create randomize run for each individual neuron of the ensemple
%--------------------------------------------------------------------------
%D0:
%-Manually enter: powerVector, durationVector
%--------------------------------------------------------------------------

%TODO: streamline the saving of stim parameters
%% Runs connectivity
%define where to save the file
if pl.Connected()
    pl.Disconnect();
end
savePrairieFiles(savePath, pl, 'connectivity_pre')

%Manually Enter: 
powerVector = .004*[15 15 35];
durationVector = [20 20 20];

% powerVector = .004*[15 50 15 15];
% durationVector = [20 20 20 20];

spiralVector = [20 20 20];
% spiralVector = [20 20 20 20];
initDelay = 5000;
reps = 5;
numberNeurons=3; %4
iterations = 1;

createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector,iterations, 'connectivity_pre_', true)
% loadCommand = "-lmp " + fullfile(savepath, "connectivty_holostim.xml");

loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);
%% Run 
%--------------------------------------------------------------------------
%D0:
%1) (zero PMT, power) Add water to imaging window if needed
%2) Load xml
%--------------------------------------------------------------------------
clear s
ConnectivityAcqnvsPrairie(folder, animal, day, StimMask, 'PRE')


%% create stims for pretrain
initDelay = 0;
reps = 1;
numberNeurons=3;%4
iterations = 121;
createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector, iterations, 'preTrain', false)

%%
%Compute vectorHolo
%--------------------------------------------------------------------------
%D0:
%1) Confirm IHSI mean, range
%2) seedBase - if we will seed the baseline, then set to 1.  
% - if seedBase 0, we wait for baseline before starting stims
%--------------------------------------------------------------------------

frameRate = 30 %baseline_frameRate
baseFrames = 2*60*30; 
expectedLengthExperiment = 40*60*frameRate

% IHSImean, IHSIrange
IHSImean = 20; 
IHSIrange = 10; 
[vectorHolo, ISI] = createVectorHolo(frameRate, expectedLengthExperiment, IHSImean, IHSIrange, false);

seedBase = 0; %Set this to 1 if you will seed the baseline
if ~seedBase
    vectorHolo = vectorHolo + baseFrames;
end
% num imaging reps should be 75600 = 72000+3600

%% create masks bot and image to check during experiment
createBot(savePath, x(ensembleNeurons),y(ensembleNeurons))
ensembleMask = holoMaskRedGreen;
ensembleMask(~ismember(ensembleMask,ensembleNeurons))= 0;
figure('Position', [600,300, 256, 256])
imshow(ensembleMask);

%% run Pre-training

%%
%Seed BMI baseVal using Pretrain
%--------------------------------------------------------------------------
%D0:
%1) seedBase - if we will seed the baseline, then set to 1. 
% - if seedBase 0, we wait for baseline before starting stims
%2) Copy-paste BMI_target_info filename (into 'pretrain_file')
%--------------------------------------------------------------------------
seedBase = 0; 
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
if seedBase
    %TODO:
    pretrain_file = 'BMI_online190523T010653'
    load(fullfile(savePath, pretrain_file)); 
    pretrain_base = data.baseVector; 
    pretrain_base(:, isnan(pretrain_base(1,:))) = [];
    baseValSeed = pretrain_base(:,end)
end

%%
% Pre-training

%--------------------------------------------------------------------------
%D0:
%
%Copy Paste BMI_target_info_.mat
%Change the Mark Points:
%Clear Point Series, Load pretrain.xml
%Make 121 iterations, put "Wait for Trigger" = First Reptition, Trigger
%Trigger Selection: Start with External, PFI1
%
%Change BOT to take only the E2 ensemble.  This lets us see the E2 neurons
%better, online.

% Then, before running cell:
%1) start video
%2) start load cells
%3) start pyctrl
%--------------------------------------------------------------------------

clear s
baselineCalibrationFile = 'BMI_target_info_20190526T113926.mat';
vectorVTA = []
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

expt_str = 'HoloVTA_pretrain'; 
debug_bool = 0; 
debug_input = []; 


BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input, baseValSeed);
%saves filename with expt_str
% BMIAcqnvsPrairienoTrialsHoloCL(folder, animal, day, expt_str, baselineCalibrationFile, baseline_frameRate, vectorHolo, vectorVTA, cursor_zscore_bool)


%35 holo hits from pretrain 1
%

%--------------------------------------------------------------------------
%D0:
%Stop:
%1) pyctrl
%2) load cells
%3) video
%--------------------------------------------------------------------------

%% run BMI
%--------------------------------------------------------------------------
%D0:
%remember to set the markpoints to proper stim. Remove mark points!!!!
% Then, before running cell:
%0) put water under objective
%1) start video
%2) start load cells
%3) start pyctrl
%--------------------------------------------------------------------------


%Get baseValSeed from HoloVTA_pretrain!  load file, take the last valid
%baseVal
% load_baseVal = 0; 
% if load_baseVal
% baselineCalibrationFile = 'BMI_target_info_20190523T220638.mat';
pretrain_file = 'BMI_online190524T131817'
load(fullfile(savePath, pretrain_file)); 
pretrain_base = data.baseVector; 
pretrain_base(:, isnan(pretrain_base(1,:))) = [];
baseValSeed = pretrain_base(:,end)

%%
% baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
% baselineCalibrationFile = 'BMI_target_info_20190523T110626.mat';
%%
%  baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
vectorHolo = [];
vectorVTA= []; 
debug_bool = 0; 
debug_input = []; 
cursor_zscore_bool = 0; 
BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, ...
    'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input, baseValSeed)

%Stop:
%1) pyctrl
%2) load cells
%3) video


%% Holo stim checking connectivity
if pl.Connected()
    pl.Disconnect();
end
% create randomize run for each individual neuron of the ensemple
savePrairieFiles(savePath, pl, 'connectivity_post')
% numberNeurons, reps, power, varName, flagRandom)

powerVector = .004*[20 50 28 15];
durationVector = [20 40 20 20];
spiralVector = [20 20 20 20];
initDelay = 5000;
reps = 5;
numberNeurons=4;
iterations = 1;

createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector, iterations, 'connectivity_post_', true)

pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_all.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%% run

clear s
ConnectivityAcqnvsPrairie(folder, animal, day, ensembleMask, 'POST')
%% end


%%
%--------------------------------------------------------------------------
%D0:

%1) SAVE THE WORKSPACE IN FOLDER
%2) SAVE THIS Protocol script IN FOLDER (savePath)
%3) Start converting imaging data
%3) Remove mouse
%4) collect load cell sensor baseline data

%--------------------------------------------------------------------------
%%
%NOTES:
%first test of E1 threshold being more strict
%did pretrain and fucked the cells up, bmi was weak can analyze effect of
%fucked up cells