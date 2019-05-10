%%%TODO'S!!!!
%{
- Create environments for each part of the experiment
- holostim trigger
- frame trigger
%}

%% Main protocol for the experiment


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

%% Select ChroME neurons
% Select area of interest with the 2p
% find red neurons
[holoMask, Im, Img, px, py] = makeMasksPrairie();

% it returns the mask that will be used for the holostim

% if we want to add neurons
[holoMask, ~,~] = addcell (Im, holoMask,9);
% delete neurons that we don't want by position on image
holoMask = deleteMask(Im, holoMask, 3);  % third var is the number of areas to delete

%% when we are happy with the result
[x,y] = findCenter(holoMask);
red = [x';y'];
Im = double(Im);
Img = double(Img);

% save red and Im in folder/animal/day
filetosave = fullfile(savePath, 'red.mat');
save(filetosave,'Im', 'Img', 'red', 'holoMask')

%% prepare HOLO STIM 
% MAKE SURE YOU DO NOT SAVE RED CHANNEL HERE!!! 
% load environment
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo_all.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

% creates holos with the mask of the red components as input
createGplFile(savePath, holoMask, posz, px)
% upload the .gpl file in the SLM 

%define where to save the file
savePrairieFilesHolo(savePath)

% create the stim train for prairie
createXmlFile(savePath, max(max(holoMask)), 1, false)

%% Run HOLO STIM

HoloAcqnvsPrairie(folder, animal, day, holoMask)


%% Obtain spatial components
% run OnAcidPrairieNotebook.ipynb 

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
% load the VoltageRec to check the results of holoStim
plotHoloStimTimeLock(holoActivity, voltageRec, wd) % --> To plot the result of

%% Baseline acquisition

% loads the result of OnAcid
load(fullfile(savePath,'redcomp.mat'));

% runs the baseline acquisiton
% reminder to remove the RED channel
% Baseline environment already removes MARKPOINTS and set the reps to 27000
BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMask, onacid_bool);
% TODO SYNCHRONIZATION WITH PRAIRIE
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

% load by hand! --> (you can blame Vivek for this :P load(fullfile(savePath,'BaselineOnline.mat'));

%% Selection of neurons
% plots neurons so we can select which ones we like the most 

if onacid_bool
    totalneurons = min(size(AComp,2), 20);
else
    totalneurons = max(max(holoMask));
end

plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)

%% Baseline simulation
% select correct parameters on
vivek_tb_test_baseline_to_calibration

% run the simulation of baseline
% baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, ...
%     target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 

%% Holo stim checking 4 together
% TODO make this smoother, by being able to plot results + creating the
% file + creating the xml etc
% create the masks for holo stim the 4 ensemble neurons
if onacid_bool
    EnsembleMask=deleteNonEnsemble (AComp, E2_base, px, py);
else
    EnsembleMask = zeros(size(holoMask,1));
    for indn = 1:length(E2_base)
        auxmask = holoMask;
        auxmask(auxmask~=E2_base(indn)) = 0;
        auxmask(auxmask~=0) = indn;
        EnsembleMask = auxmask + EnsembleMask;
    end
end
savePrairieFiles(savePath, pl, '4neurons_together')
createGplFile(savePath, EnsembleMask, posz, px, 'ensemble_')

% upload the GPL file
% run 4 masks together
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%% Holo stim checking connectivity
% create randomize run for each individual neuron of the ensemple
savePrairieFiles(savePath, pl, 'connectivity_pre')
createXmlFile(savePath, 4, 5, 'connectivty_', true)
% ran each neuron independently

%%
%Compute vectorHolo
frameRate = 30 %baseline_frameRate
expectedLengthExperiment = 40*60*frameRate
vectorHolo = createVectorHolo(baseline_frameRate, expectedLengthExperiment);

%%
%run Pre-training
clear s
baselineCalibrationFile = 'BMI_target_info.mat';
vectorVTA = []
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

expt_str = 'HoloVTA_pretrain'; 
BMIAcqnvsPrairienoTrials(folder, animal, day, expt_str, baselineCalibrationFile, baseline_frameRate, vectorHolo, vectorVTA)

%% run BMI
%remember to set the markpoints to proper stim
BMIAcqnvsPrairienoTrials(folder, animal, day, 'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)

%% Holo stim for functional connectivity
% upload the GPL file
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile('F:/VivekNuria/utils', "Tseries_VivekNuria_holo.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%% Holo stim checking connectivity
% create randomize run for each individual neuron of the ensemple
savePrairieFiles(savePath, pl, 'connectivity_pre')
createXmlFile(savePath, 4, 5, 'connectivty_', true)
% ran each neuron independently
 