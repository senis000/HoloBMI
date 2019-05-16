%%%TODO'S!!!!
%{
%put all the plots in a folder of plots
%}

%% Main protocol for the experiment
%Start by running:
%test_bmi_nidaq_triggers.m

% define Animal, day and folder where to save
animal = 'NY27'; day = 'D1';
folder = 'F:/VivekNuria/expt/HoloBmi';

% define posz TODO can we get this from prairie?
posz = 0;
frameRate = 29.989;

onacid_bool = false;

savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%%
%Once the imaging view is set, disable the motor control!!!!

%% Select ChroME neurons
% Select area of interest with the 2p
% find red neurons
[holoMask, Im, Img, px, py] = makeMasksPrairie();

% it returns the mask that will be used for the holostim

% if we want to add neurons
[holoMask, ~,~] = addcell (Im, holoMask,9); %Press enter in beginning


holoMask = deleteMask(Im, holoMask, 3);  % third var is the number of areas to delete
%THIS DOESNT NEED ENTER AT BEGINNING
% delete neurons that we don't want by position on image

% if we want to add neurons in green channel
[holoMaskRedGreen, ~,~] = addcell (Img, holoMask,9);
%DONT DELETE NEURONS HERE!  THEY ARE FROM RED CHANNEL


%% when we are happy with the result
[x,y] = findCenter(holoMask);
[xrg,yrg] = findCenter(holoMaskRedGreen);
red = [x';y'];
redGreen = [xrg';yrg'];
Im = double(Im);
Img = double(Img);

% save red and Im in folder/animal/day
filetosave = fullfile(savePath, 'red.mat');
save(filetosave,'Im', 'Img', 'red', 'redGreen', 'holoMask', 'holoMaskRedGreen')

%% prepare HOLO STIM 
% MAKE SURE YOU DO NOT SAVE RED CHANNEL HERE!!! 
% Change number of reps depending on the amount of neurons in prairie view

% load environment
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_all.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

% creates holos with the mask of the red components as input
createGplFile(savePath, holoMask, posz, px)
% upload the .gpl file in the SLM 
%(Import: Top half of MarkPoints)

%define where to save the file
savePrairieFilesHolo(savePath)

% create the stim train for prairie
createXmlFile(savePath, max(max(holoMask)), 1, 0.3, '', false)
%Import: Mark Point Series, 

%% Update prairie view repetitions based on num neurons to stim
num_stim_neurons = max(max(holoMask));
stim_time_per_neuron = 2.1;
disp(num_stim_neurons*stim_time_per_neuron*frameRate)

%% Run HOLO STIM
%This does one neuron at a time.
clear s
HoloAcqnvsPrairie(folder, animal, day, holoMask)
%WE WANT RED+GREEN FOR THIS

%TODO: make this closed loop, and wait for the neurons to be inactive
%before stimming them.

%%
%Convert the holostim acqn using image-block ripping utility
%Load holoactivity.mat and voltage recording for plotting
%(Import to matlab: output type is Numeric Matrix.  Name it "voltageRec")
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

%% Baseline acquisition
% loads the result of OnAcid
% REMOVE RED CHANNEL!!!!

%First:
%1) start video
%2) start load cells
%3) start pyctrl

if ~onacid_bool
    AComp = 0;
else
    load(fullfile(savePath,'redcomp.mat'));
end

% runs the baseline acquisiton
% reminder to remove the RED channel
% Baseline environment already removes MARKPOINTS and set the reps to 27000
BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMaskRedGreen, onacid_bool, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

%Stop:
%1) pyctrl
%2) load cells
%3) video

%% Selection of neurons
% plots neurons so we can select which ones we like the most 

% load by hand! --> (you can blame Vivek for this :P load(fullfile(savePath,'BaselineOnline.mat'));
if onacid_bool
    totalneurons = min(size(AComp,2), 20);
else
    totalneurons = max(max(holoMaskRedGreen));
    CComp = [];
    YrA = []; 
end

%Copy paste base_file path: 
base_file = fullfile(savePath, 'BaselineOnline190515T214233.mat')
load(base_file); 
plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)

%%
%Manually enter:
E1_base = sort([7 13 10 4], 'ascend') %JUST NEEDS GCAMP
E2_base = sort([1 27 25  10], 'ascend') %NEEDS CHROME

% E2_candidates = [39 45 59 37 88 6 26 46 78 48 22 20 33]
%TODO: 
%onacid from the baseline acquisition?
%if we don't do onacid on the baseline activity, we won't have any
%roi's just from gcamp

%%
%If you need to use BMI data as the baseline data: 
% bmi_file = fullfile(savePath, 'BMI_online190515T010526.mat'); 
% bmi_data = load(bmi_file); 
% bmi_base = fullfile(savePath, ['base_' 'BMI_online190515T010526.mat']);
% baseActivity = bmi_data.data.bmiAct(:, ~isnan(bmi_data.data.bmiAct(1,:))); 
% save(bmi_base, 'baseActivity'); 
% 
% E1_base = [1 2 3 4]; 
% E2_base = [5 6 7 8]; 

%% Calibrate Target with Baseline simulation
% select correct parameters on
% vivek_tb_test_baseline_to_calibration


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

sec_per_reward_range = [100 80]; 

frames_per_reward_range = sec_per_reward_range*baseline_frameRate %[1 1.5]*60*frameRate
%multiply by frames per minute to convert
%to 

target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 4*60*ceil(frameRate)
dff_win_bool = 1
dff_win = 2
 
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

% baseline2target(n_f_file, A_file, onacid_bool, E1_base, E2_base, frames_per_reward_range, ...
%     target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, savePath)


% run the simulation of baseline
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 

%%
%copy paste the path to the target info file
target_info = load(fullfile(savePath, 'BMI_target_info.mat'));

%% Holo stim checking 4 together
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
pl = actxserver('PrairieLink.Application');
pl.Connect();
savePrairieFiles(savePath, pl, '4neurons_together')
createGplFile(savePath, EnsembleMask, posz, px, 'ensemble_')

createBot(savePath, x(E2_base),y(E2_base))

loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);

%ACTIONS:
% SELECT BOT on image window, import from file
% upload the GPL file
% Manually sselect neurons, make group, laser power 50, 4 iterations with
% 10 sec bw
% CHANGE NUMBER OF REPS IN T-SERIES
%% Run 4 masks together

pl.SendScriptCommands("-ts"); 
pause(1);

%% Holo stim checking connectivity
% create randomize run for each individual neuron of the ensemple
savePrairieFiles(savePath, pl, 'connectivity_pre')
power = 0.1; %0.2 = 50
createXmlFile(savePath, 4, 5, power, 'connectivity_pre_', true)
% loadCommand = "-lmp " + fullfile(savepath, "connectivty_holostim.xml");
% pl.SendScriptCommands(loadCommand);
% LOAD XML FILE

%% Run 
pl.SendScriptCommands("-ts"); 
pause(1);

%% disconnect prairie
pl.Disconnect()
% ran each neuron independently

%%
%Compute vectorHolo
frameRate = 30 %baseline_frameRate
expectedLengthExperiment = 40*60*frameRate
[vectorHolo, ISI] = createVectorHolo(frameRate, expectedLengthExperiment, 20, 5, false);

%TODO: give input: baseFrames, vectorHolo will only generate stims after
%baseFrames

%%
%run Pre-training
%Change the Mark Points:
%Make 120 reps, put "Wait for Trigger" = First Reptition, Trigger
%Selection: Start with External, PFI1

% Then, before running cell:
%1) start video
%2) start load cells
%3) start pyctrl

clear s
baselineCalibrationFile = 'BMI_target_info_20190515T214817.mat';
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

baseValSeed = ones(length(E1_base)+length(E2_base), 1)+nan
BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input, baseValSeed)
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
pretrain_file = ''
load(fullfile(savePath, pretrain_file)); 
data.baseVal

pretrain_base = data.baseVector; 
pretrain_base(isnan(pretrain_base(1,:))) = [];

BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, ...
    'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input)

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
createXmlFile(savePath, 4, 5, 0.2, 'connectivity_post_', true)% ran each neuron independently

%% run
pl.SendScriptCommands("-ts"); 
pause(1);
%% end
pl.Disconnect()

%%
%SAVE THIS IN FOLDER

 