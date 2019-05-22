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
animal = 'NY28'; day = 'D3_p2';
folder = 'E:/VivekNuria/expt/HoloBmi';

% define posz TODO can we get this from prairie?
posz = 0;
frameRate = 29.989;

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
%--------------------------------------------------------------------------
%% Select ChroME neurons
% Select area of interest with the 2p
% find red neurons

%--------------------------------------------------------------------------
%DO:
%-mainMasksPrairie, make holoMask
%-deleteMask, delete bad neurons.  
%-addcell, add neurons
%-add green cells: 
%--[holoMaskRedGreen, ~,~] = addcell (Img, holoMask,9);
%--[holoMaskRedGreen, ~,~] = addcell (Img, holoMaskRedGreen,9);
%--------------------------------------------------------------------------


[holoMask, Im, Img, px, py] = makeMasksPrairie();

% it returns the mask that will be used for the holostim

% if we want to add neurons
[holoMask, ~,~] = addcell (Im, holoMask,9); %Press enter in beginning


holoMask = deleteMask(Im, holoMask, 7);  % third var is the number of areas to delete
%THIS DOESNT NEED ENTER AT BEGINNING
% delete neurons that we don't want by position on image

% if we want to add neurons in green channel
[holoMaskRedGreen, ~,~] = addcell (Img, holoMask,9);
%If you need to do more than one round:
[holoMaskRedGreen, ~,~] = addcell (Img, holoMaskRedGreen,9);
%DONT DELETE NEURONS HERE!  THEY ARE FROM RED CHANNEL

%%
%See the holoMask:
h = figure;
imshow(holoMask)

%% Save the holoMask
[x,y] = findCenter(holoMask);
[xrg,yrg] = findCenter(holoMaskRedGreen);
red = [x';y'];
redGreen = [xrg';yrg'];
Im = double(Im);
Img = double(Img);

% save red and Im in folder/animal/day
filetosave = fullfile(savePath, 'red.mat');
save(filetosave,'Im', 'Img', 'red', 'redGreen', 'holoMask', 'holoMaskRedGreen')

%% prepare HOLO STIM of individual neurons
%Create gpl, xml files for individual points

% load environment
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_all.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

% creates holos with the mask of the red components as input
createGplFile(savePath, holoMask, posz, px)


%define where to save the file
savePrairieFilesHolo(savePath)

% create the stim train for prairie
num_neurons_stim = max(max(holoMask));

numberNeurons=max(max(holoMask));
power = 0.2; 
duration = 30; 
numSpirals = 20; 

powerVector = power*ones(1,numberNeurons); %*0.1; %0.2 = 50, 0.1=25
durationVector = duration*ones(1,numberNeurons);
spiralVector = numSpirals*ones(1,numberNeurons);
initDelay = 2000;
iterations = 1;
reps = 1;

createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector, iterations, '', false)

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

clear s
HoloAcqnvsPrairie(folder, animal, day, holoMask)
%WE WANT RED+GREEN FOR THIS
%TODO: make this closed loop, and wait for the neurons to be inactive
%before stimming them.

%--------------------------------------------------------------------------
%DO: 
%Image-block ripping utility: Convert the holostim acqn
%Load holoactivity.mat and voltage recording for plotting
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
plotHoloStimTimeLock(holoActivity, voltageRec, 40, 100) % --> To plot the result of
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
% E2_base = sort([158 74 77 147 113 149 65 191 60 108 69 98 116 164 161 129 26], 'ascend')

% E2_base =[26    60    65    69    74    77 ...
%     98   108   113   116   129   147 ...
%    149   158   161   164]; 

E2_base = [98   158   161   164]; 
% E2_base = sort([14 10 26 9], 'ascend')

%TODO: update protocol to iterate from E2_candidates to E2_base

%% Holo stim of Ensemble
%TODO create file + creating the xml etc
% create the masks for holo stim the 4 ensemble neurons
if onacid_bool
    EnsembleMask=deleteNonEnsemble (AComp, E2_base, px, py);
else
    EnsembleMask = zeros(size(holoMaskRedGreen,1));
    for indn = 1:length(E2_base)
        auxmask = holoMaskRedGreen;
        auxmask(auxmask~=E2_base(indn)) = 0;
        auxmask(auxmask~=0) = indn;
        EnsembleMask = auxmask + EnsembleMask;
    end
end
figure;
imshow(EnsembleMask); 

%%
%GPL for Stim Ensemble
pl = actxserver('PrairieLink.Application');
pl.Connect();
savePrairieFiles(savePath, pl, '4neurons_together')
createGplFile(savePath, EnsembleMask, posz, px, 'ensemble_')

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
%1
% p15,d20
% pretty good, missed a couple stims, afraid it will green out
%2 
% p15, 220
% very good, but very active
%
%3
%ABORT
%
%4
%p40,d20
%pretty good, worried about over green
%
%5
%ABORT
%
%6
%ABORT
%
%7
%Good
%p15,d20
%
%8
%NO
%
%9
%NO just weird
%
%10
%NO
%
%11
%NO
%
%12
%it seems to respond but veryyyy slowly...
%
%13
%NO
%
%14
%VERY NICE
%p22, d20
%
%15
%VERY NICE
%p50, d30
%
%16
%GOOD
%p40, d30

powerVector = [15 22 50 40];
durationVector = [20 20 30 30];

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
%0) put water
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
BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMaskRedGreen, onacid_bool, frameRate);
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
base_file = fullfile(savePath, 'BaselineOnline190521T222151.mat')
load(base_file); 
% plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)
plotNeuronsBaseline(baseActivity, CComp, YrA, 20)

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
E1_base = sort([64 19 141 83], 'ascend') %JUST NEEDS GCAMP
% E2_base = sort([7 12 34 48], 'ascend')


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

sec_per_reward_range = [120 100]; 

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

 baseline2target_vBMI(n_f_file, A_file, onacid_bool,  ...
    E1_base, E2_base, frames_per_reward_range, target_on_cov_bool, ...
    prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath, ...
    cursor_zscore_bool, f0_init_slide);
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
% T: 0.6
% num_valid_hits: 8
% num_hits: 102
%--------------------------------------------------------------------------


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

savePrairieFiles(savePath, pl, 'connectivity_pre')
% powerVector = [0.06 0.04 0.08 0.06]; %0.2 = 50

%Manually Enter: 
powerVector = .004*[15 22 50 40];
durationVector = [20 20 30 30];
spiralVector = [20 20 20 20];
initDelay = 5000;
reps = 5;
numberNeurons=4;
iterations = 1;

createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector,iterations, 'connectivity_pre_', true)
% loadCommand = "-lmp " + fullfile(savepath, "connectivty_holostim.xml");

loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);
%% Run 
%--------------------------------------------------------------------------
%D0:
%1) Add water to imaging window if needed
%--------------------------------------------------------------------------

pl.SendScriptCommands("-ts"); 
pause(1);

%% disconnect prairie
pl.Disconnect()

%% create stims for pretrain
initDelay = 0;
reps = 1;
numberNeurons=4;
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

seedBase = 1; %Set this to 1 if you will seed the baseline
if ~seedBase
    vectorHolo = vectorHolo + baseFrames;
end
% num imaging reps should be 75600 = 72000+3600

%%
%Seed BMI baseVal using Pretrain
%--------------------------------------------------------------------------
%D0:
%1) seedBase - if we will seed the baseline, then set to 1. 
% - if seedBase 0, we wait for baseline before starting stims
%2) Copy-paste BMI_target_info filename (into 'pretrain_file')
%--------------------------------------------------------------------------
seedBase = 1; 
baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan;
if seedBase
    %TODO:
    pretrain_file = 'BMI_online190521T224159'
    load(fullfile(savePath, pretrain_file)); 
    pretrain_base = data.baseVector; 
    pretrain_base(:, isnan(pretrain_base(1,:))) = [];
    baseValSeed = pretrain_base(:,end)
end

%%
% Pre-training

%Copy Paste BMI_target_info_.mat
%Change the Mark Points:
%Clear Point Series, Load pretrain.xml
%Make 120 reps, put "Wait for Trigger" = First Reptition, Trigger
%Selection: Start with External, PFI1

% Then, before running cell:
%1) start video
%2) start load cells
%3) start pyctrl

clear s
baselineCalibrationFile = 'BMI_target_info_20190521T224235.mat';
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


%Stop:
%1) pyctrl
%2) load cells
%3) video


%%

%% run BMI
%remember to set the markpoints to proper stim. Remove mark points!!!!

% Then, before running cell:
%1) start video
%2) start load cells
%3) start pyctrl

%Get baseValSeed from HoloVTA_pretrain!  load file, take the last valid
%baseVal
% load_baseVal = 0; 
% if load_baseVal
baselineCalibrationFile = 'BMI_target_info_20190521T224235.mat';
pretrain_file = 'BMI_online190521T232745'
load(fullfile(savePath, pretrain_file)); 
pretrain_base = data.baseVector; 
pretrain_base(:, isnan(pretrain_base(1,:))) = [];
baseValSeed = pretrain_base(:,end)

%%
% baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
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
%% Holo stim for functional connectivity
% upload the GPL file
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_all.env");
pl.SendScriptCommands(loadCommand);

%% Holo stim checking connectivity
% create randomize run for each individual neuron of the ensemple
savePrairieFiles(savePath, pl, 'connectivity_post')
% numberNeurons, reps, power, varName, flagRandom)

powerVector = .004*[15 22 50 40];
durationVector = [20 20 30 30];
spiralVector = [20 20 20 20];
initDelay = 5000;
reps = 5;
numberNeurons=4;
iterations = 1;

createXmlFile(savePath, numberNeurons, reps, initDelay, durationVector, powerVector, spiralVector, iterations, 'connectivity_post_', true)

%% run
pl.SendScriptCommands("-ts"); 
pause(1);
%% end
pl.Disconnect()

%%
%SAVE THIS Protocol script IN FOLDER (savePath)
%collect load cell sensor baseline data

%%
%NOTES:
%In holostim video:
%between 8500 and 9500 frames, I saw a wave of activity, maybe a seizure?
%This video is interesting to analyze.
%TODO: We need this stim to be closed loop so we don't waste stims...
%
%CHECK cell14 video
 %
%One idea is to deliver low intensity stims until the target pattern is
%reached, so we ensure a stim eventually leads to target hit.
%
%We should add time between rewards so one large transient can't cause
%multiple rewards?

%Noticed a motion before 26900 frames, so i adjusted manually.  
%I'm running 40000 frames on F.  I'll have to run 32000 frames onto E

%
%can i script the image-block ripping utility?  
%can i script the loadcell conversion?
%%
%BMIp2
%BMIp1
%pretrain

%frame 25590, in BMIp2, the animal isn't even licking the sucrose.  it's
%hanging on the spout.
