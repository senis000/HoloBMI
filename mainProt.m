%%%TODO'S!!!!
%{
- We need to know which voltage-output is which
- plotHoloStimTimeLockFromTiff --> adjust voltage-output index and test
- Make sure we know how to obtain the BoT with holo+voltage
- Create environments for each part of the experiment
- adapt code to match neurons (INDEX) from onacid with neurons from the holomask
- Synchronization with Prairie -> each frame in voltage-ouput
- check port, pin and pause time for arduino
- holostim trigger
- frame trigger
%}

%% Main protocol for the experiment


% define Animal, day and folder where to save
animal = 'NY27'; day = 'D1';
folder = 'F:/VivekNuria/expt/HoloBmi';

% define posz TODO can we get this from prairie?
posz = 176.5;
frameRate = 29.989;

savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%% Select ChroME neurons
% Select area of interest with the 2p
% find red neurons
[holoMask, Im, Img, px, py] = makeMasksPrairie();

% it returns the mask that will be used for the holostim

% if we want to add neurons
[holoMask, ~,~] = addcell (Im,holoMask,9);
% delete neurons that we don't want by position on image
holoMask = deleteMask(Im, holoMask, 3);  % third var is the number of areas to delete

%% when we are happy with the result
[x,y] = findCenter(holoMask);
red = [x';y'];
Im = double(Im);

% save red and Im in folder/animal/day
filetosave = fullfile(savePath, 'red.mat');
save(filetosave,'Im', 'Img', 'red', 'holoMask')

%% prepare HOLO STIM 
% MAKE SURE YOU DO NOT SAVE RED CHANNEL HERE!!! 
% creates holos with the mask of the red components as input
createGplFile(savePath, holoMask, posz, px)
% creates regions of interest with the mask of the red components as input
createBot(savePath, x,y)
% upload the .gpl file in the SLM and the botfile in the BoT

%define where to save the file
savePrairieFilesHolo(savePath)
%%
createXmlFile(savePath, max(max(holoMask)), 1, false)

%% Run HOLO STIM
% TODO set the environment for the holo recording
%check -TSeriesLoad (-tsl) ["path"] with a saved environment

% load the files BoT and VoltageRec to check the results of holoStim
% plotHoloStimTimeLock(botData, voltageRec, wd)  --> To plot the result of
% stim

%% Obtain spatial components
% run OnAcidPrairieNotebook.ipynb 
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
% load the files BoT and VoltageRec to check the results of holoStim
plotHoloStimTimeLock(botData, voltageRec, wd)


%% Baseline acquisition

% loads the result of OnAcid
load(fullfile(savePath,'redcomp.mat'));

% runs the baseline acquisiton
BaselineAcqnvsPrairie(folder, animal, day, AComp);
% TODO SYNCHRONIZATION WITH PRAIRIE
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

load(fullfile(savePath,'BaselineOnline.mat'));

%% Selection of neurons
% plots neurons so we can select which ones we like the most 
totalneurons = min(size(AComp,2), 20);
plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)

%% Baseline simulation
% select correct parameters on
vivek_tb_test_baseline_to_calibration

% run the simulation of baseline
% baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, ...
%     target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 

%% Holo stim for functional connectivity
% create the masks for holo stim the 4 ensemble neurons
EnsembleMask=deleteNonEnsemble (AComp, E2_base, px, py);
createGplFile(savePath, EnsembleMask, posz, px, 'ensemble_')
% upload the GPL file
% run 4 masks together
% create randomize run for each individual neuron of the ensemple
createXmlFile(savePath, numberNeurons, reps, flagRandom)
% upload the XML file
% ran each neuron independently

%% run Pre-training
baselineCalibrationFile = 'BMI_target_info.mat';

BMIAcqnvsPrairienoTrials(folder, animal, day, expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)

%% run BMI
BMIAcqnvsPrairienoTrials(folder, animal, day, 'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)

%% Holo stim for functional connectivity
% ran each neuron independently
 