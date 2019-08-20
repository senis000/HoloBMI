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

% define Animal, day and folder where to save
animal = 'NY107'; day = 'D1';
folder = 'E:\Vivek_e\training';
% folder = 'F:\Vivek\training';

% define posz TODO can we get this from prairie?
posz = 0;
savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%Confirm the BMI settings in this function: 
[task_settings] = define_BMI_task_settings(); 
if(task_settings.fb.fb_bool) 
    a = arduino(task_settings.fb.arduino.com, ...
        task_settings.fb.arduino.label);
else
    a = [];
end

screen_size = get(0,'ScreenSize');

%%
%--------------------------------------------------------------------------
%DO:
%find FOV
%-disable the motor control!!!!
%-set imaging FOV.  
%-ACTIVATE RED+GREEN channels
%-collect 1000 frame video
%--------------------------------------------------------------------------
%DO:
%-convert with Image-Block Ripping Utility
%-drag into ImageJ
%-split into red and green channels
%--Image>Stacks>Tools>Make Substack
%-take average frame of the two channels separately:
%--Image>Stacks>Z Project 
%--save the images into 'redgreen' folder in 'savePath'
%--------------------------------------------------------------------------
%%
%Load red, green images: 
redgreen_dir = fullfile(savePath, 'redgreen'); 
exist(redgreen_dir)
red_path      = fullfile(redgreen_dir, 'red_mean.tif'); 
exist(red_path)
green_path    = fullfile(redgreen_dir, 'green_mean.tif'); 
exist(green_path)
%%
% red_path      = fullfile('E:\ines_e\redgreen-000', 'red.tif'); 
% green_path    = fullfile('E:\ines_e\redgreen-000', 'green.tif'); 

green_im    = imread(green_path); 
red_im      = imread(red_path);  

%Initialize data structure to save overlays:
rg_struct = struct(...
    'im', [], ...
    'rg_minmax', [], ...
    'rg_minmax_perc', [], ...
    'r_min', [], ...
    'r_min_perc', [], ...
    'g_min', [], ...
    'g_min_perc', [], ...
    'r_max', [], ...
    'r_max_perc', [], ...
    'g_max', [], ...
    'g_max_perc', []...
    ); 
num_overlays = 0;

%%
[rg_struct, num_overlays] = pick_rg_overlay(red_im, green_im, rg_struct, num_overlays);

%%
%--------------------------------------------------------------------------
%DO:
%(Choose an overlay image for drawing ROIs.)
%Select i
%--------------------------------------------------------------------------
i = num_overlays; 
im_bg = rg_struct(i).im;
h = figure; imagesc(im_bg); axis square; 

%% Select Red +  Green neurons
plot_images = struct('im', [], 'label', ''); 
plot_images(1).im = green_im; 
plot_images(1).label = 'Green Mean'; 
plot_images(2).im = red_im; 
plot_images(2).label = 'Red Mean'; 

roi_init_bool = 1; %Initializes an empty roi_data structure
if(roi_init_bool)
    roi_data = init_roi_data(im_bg); 
end
disp('Adding ROIs to image!'); 
[roi_data] = draw_roi_2chan(plot_images, roi_data);

%%
%Add more ROI if needed: 
[roi_data] = draw_roi_2chan(plot_images, roi_data);
%%
%Delete ROI if needed: 
[roi_data] = delete_roi_2chan(plot_images, roi_data);
%%
%See the roi_mask:
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(roi_data.im_roi_rg); axis square; title('roi mask'); 

h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
imagesc(roi_data.roi_mask); axis square; title('roi mask');    


%% Compute strcMask for roi_data: 



%% Save the roi_data
%plot_images
%rg_struct
%roi_data
%roi_ctr - this is needed for drawing BOT so we can see the neurons during
%the experiment. 
%
% 
%HERE
[x,y]           = findCenter(roi_data.roi_mask);
roi_ctr.xy      = [x';y']; %size: 2 x num_roi
roi_ctr.chan    = repmat(struct('xy', []), [2 1]); 
for chan_i = 1:2
%     [x,y] = findCenter(roi_data.chan(chan_i).roi_mask); 
    roi_ctr.chan(chan_i).xy = roi_ctr.xy(:,roi_data.chan(chan_i).idxs); 
end
roi_mask = roi_data.roi_mask;
%SAVE: 
roi_data_file = fullfile(savePath, 'roi_data.mat');
save(roi_data_file,'roi_mask', 'plot_images', 'rg_struct', 'roi_data', 'roi_ctr')

% save(filetosave,'Im', 'Img', 'red', 'redGreen', 'holoMask', 'holoMaskRedGreen')

% %%
% strcMask = obtainStrcMaskfromMask(roi_data.roi_mask)
% 
% %%
% strcMask_red = obtainStrcMaskfromMask(roi_data.chan(1).roi_mask)

%% create masks bot and image to check during experiment
x_roi = roi_ctr.xy(1,:); 
y_roi = roi_ctr.xy(2,:); 
bot_base_path = fullfile(savePath, 'Bot_base.cfg'); 
createBot(bot_base_path, x_roi,y_roi)

%%
% roi_ind = unique(roi_data.chan(1).roi_mask(:));
% roi_ind(roi_ind==0) = []; %remove the 0 ind
% num_roi = length(roi_ind)

%%

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

%DO THIS:
%TEST THIS!  Would rather do this than manually select ROI's.  

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

if ~task_settings.onacid_bool
    AComp = 0;
else
    %load onacid results
    load(fullfile(savePath,'roi_data_file.mat'));
end


% Baseline environment already removes MARKPOINTS and set the reps to 27000

clear s
[baseline_mat, baseline_dat] = ...
    BaselineAcqnvsPrairie(folder, animal, day, AComp, roi_mask, task_settings);

% BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMaskRedGreen, onacid_bool, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

%TODO: 
%Confirm that BaselineAcqnvsPrairie works properly with AComp from OnAcid
%Edit so that we can see Green vs Red+Green neurons

%--------------------------------------------------------------------------
%D0:
%0) Abort T-series (cuz of voltage recording)
%1) B: pyctrl stop
%2) B: load cells stop
%3) B: video stop
%4) B: Drag load cell data to folder
%5) B: Drag video to folder
%--------------------------------------------------------------------------

%% Selection of neurons
% plots neurons so we can select which ones we like the most 
%--------------------------------------------------------------------------
%D0:
%1) Copy paste base_file name 'BaselineOnlineX.mat'
%--------------------------------------------------------------------------
% load by hand! --> (you can blame Vivek for this :P load(fullfile(savePath,'BaselineOnline.mat'));
if task_settings.onacid_bool
    totalneurons = min(size(AComp,2), 20);
else
    CComp = [];
    YrA = []; 
end
load(baseline_mat); 
%Plot Red: 
totalneurons = 30; 
sel = roi_data.chan(1).idxs; 
plotSelNeuronsBaseline(baseActivity, CComp, YrA, totalneurons, sel);
title('RED'); 
%D1
%[7 1 11 2]
%%
%Plot Green:
totalneurons = 30; 
sel = roi_data.chan(2).idxs; 
plotSelNeuronsBaseline(baseActivity, CComp, YrA, totalneurons, sel);
% title('GREEN'); 
%D2
%[42 31 40 26]
%%
%--------------------------------------------------------------------------
%D0:
%1) Choose E2_base (has to increase), E1_base (has to be suppressed)
%--------------------------------------------------------------------------
%Manually enter:
%E2 green
%E1 red
E2_base = sort([1 7 26 42], 'ascend') %Activity needs to increase 
%8 21 30 45
%R G R G
%5 6 7 8
E1_base = sort([2 11 31 40], 'ascend') %Activity needs to decrease
ensembleNeurons = [E1_base, E2_base];
ensembleID = [ones(1,length(E1_base)) 2*ones(1,length(E2_base))]; 
%Change the colors of the traces: 
plotNeuronsEnsemble(baseActivity, ensembleNeurons, ensembleID)
% E2_candidates = [39 45 59 37 88 6 26 46 78 48 22 20 33]
%%
%OPTION: Use previously collected BMI data as the baseline data: 
%
% bmi_file = fullfile(savePath, 'BMI_online190515T010526.mat'); 
% bmi_data = load(bmi_file); 
% bmi_base = fullfile(savePath, ['base_' 'BMI_online190515T010526.mat']);
% baseActivity = bmi_data.data.bmiAct(:, ~isnan(bmi_data.data.bmiAct(1,:))); 
% save(bmi_base, 'baseActivity'); 


%%
% calibration_settings_bmi1 = calibration_settings; 

%%
% calibration_settings_bmi2 = calibration_settings; 

%%
% calibration_settings_bmi3 = calibration_settings;

%%
calibration_settings_bmi4   = calibration_settings;

%% Calibrate Target with Baseline simulation
%--------------------------------------------------------------------------
%D0:
%1) Parameters: 
% - sec_per_reward_range
% - f0_win (F0: how many frames to average over)
% - dff_win (F for Dff: how many frames to average over)
%--------------------------------------------------------------------------

%7.13.19
%b2base_num_samples
%cursor_zscore_bool
%dff_win / movingAverageFrames
% baseval_win / f0_win
% b2base_val (default T/2)


% base_file = fullfile(savePath, 'BaselineOnline190514T221822.mat')
% base_file = bmi_base; 

exist(baseline_mat)
n_f_file = baseline_mat;
close all;
[calibration_settings, BMI_roi_path] = baseline2target_vE1strict_v2(n_f_file, roi_data_file, task_settings, ...
    E1_base, E2_base, savePath);
%ToDo: return the filename of 


% run the simulation of baseline
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 
%--------------------------------------------------------------------------
%D0:
%Note down: 
% - T value
% T: 0.1014
%                 num_c1: 1355
%                 num_c2: 12870
%                 num_c3: 114
%     num_hits_no_b2base: 20
%         num_valid_hits: 10
%--------------------------------------------------------------------------

%% create masks bot and image to check during experiment
x_bmi = roi_ctr.xy(1,ensembleNeurons); 
y_bmi = roi_ctr.xy(2,ensembleNeurons); 
bot_bmi_path = fullfile(savePath, 'Bot_bmi.cfg'); 
createBot(bot_bmi_path, x_bmi,y_bmi)
ensembleMask = roi_data.roi_mask;
ensembleMask(~ismember(ensembleMask,ensembleNeurons))= 0;
figure('Position', [600,300, 256, 256])
imshow(ensembleMask);
title('Mask for Ensemble Neurons'); 

%%
%Seed BMI baseVal, if you already ran BMI, and need to run again
%--------------------------------------------------------------------------
%D0:
%1) seedBase - if we will seed the baseline, then set to 1. 
% - if seedBase 0, we wait for baseline before starting stims
%2) Copy-paste BMI_target_info filename (into 'pretrain_file')
%--------------------------------------------------------------------------
seedBase = 0; 
seed_file = 'BMI_online190719T163934'
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
if seedBase
    %ENTER HERE:
    load(fullfile(savePath, seed_file)); 
    pretrain_base = data.baseVector; 
    pretrain_base(:, isnan(pretrain_base(1,:))) = [];
    baseValSeed = pretrain_base(:,end)
end

%% run BMI
%--------------------------------------------------------------------------
%D0:
% Before running cell:
%0) put water under objective
%1) start video
%2) start load cells
%3) start pyctrl
%--------------------------------------------------------------------------
vectorHolo = [];
vectorVTA= []; 
debug_bool = 0; 
debug_input = []; 
expt_str = 'BMI'; 
BMIAcqnvsPrairienoTrialsHoloCL_fb_debug_enable(folder, animal, day, ...
    expt_str, calibration_settings, task_settings, a, vectorHolo, vectorVTA, ...
    debug_bool, debug_input)
%Stop:
%1) pyctrl
%2) load cells
%3) video
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
%{
First baseline didn't work because of motion drifting and losing the ROI's
The ROIs are smaller than cortex
For future, draw roi's for baseline with BOT

TODO: 
save ensemble ROIs with rg identity
show std dev for choosing ROI
PLOT ROIMASK with Traces, check that the numbers of what we pick aren't
close to eachother
center the feedback at the mode cursor.  
%Need script to tell us how off live frames are from template image.  very
important in striatum imaging.  

First BMI was mix D1,D2
Second BMI was E1=D2, E2=D1
Third BMI was E1=D1, E2=D2
    This had calibration issue where it was very difficult to get one
    ensemble active and not the other, so not sure we can interpret this
    experiment.
%Fourth BMI was same as first.
%
%
%}