% Main protocol for the experiment

% define Animal, day and folder where to save
animal = 'a1'; day = 'd2';
folder = 'F:/VivekNuria/';

% define posz TODO can we get this from prairie?
posz = 105.25;

savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%% Select ChroME neurons
% Select area of interest with the 2p
% find red neurons
[holoMask, Im, px, py] = makeMasksPrairie(channel);

% it returns the mask that will be used for the holostim

% if we want to add neurons
[holoMask, ~,~] = addcell (Im,holoMask,9);
% delete neurons that we don't want by position on image
holoMask = deleteMask(Im, holoMask, 3);  % third var is the number of areas to delete

%when we are happy with the result
[x,y] = findCenter(holoMask);
red = [x';y'];
Im = double(Im);

% save red and Im in folder/animal/day
filetosave = fullfile(savePath, 'red.mat');
save(filetosave,'Im', 'red')

%% HOLO STIM 
% MAKE SURE YOU DO NOT SAVE RED CHANNEL HERE!!! 
% creates holos with the mask of the red components as input
createGplFile(savePath, holoMask, posz, px)
% upload the .gpl file in the SLM

%define where to save the file
%TODO DEFINE LENGTH OF EXPERIMENT DEPENDING ON NUMBER OF NEURONS
savePrairieFilesHolo(savePath)

%% Obtain spatial components

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

%% Baseline acquisition
% runs the baseline acquisiton
BaselineAcqnvsPrairie(folder, animal, day, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

%% Selection of neurons
% plots neurons so we can select which ones we like the most 
plotNeuronsBaseline(baseActivity, CComp, YrA, totalNeurons)

%% Baseline simulation
baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, ...
    target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 

%% Holo stim for functional connectivity
% create the masks for holo stim the 4 ensemble neurons
EnsembleMask=deleteNonEnsemble (AComp, E2indx, px, py);
createGplFile(savePath, EnsembleMask, posz)
% upload the GPL file
% run 4 masks together
% create randomize run for each individual neuron of the ensemple
createXmlFile(savePath, numberNeurons, reps, flagRandom)
% upload the XML file
% ran each neuron independently

%% run Pre-training
BMIAcqnvsPrairienoTrials(animal, day, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)

%% run BMI
BMIAcqnvsPrairienoTrials(animal, day, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)

%% Holo stim for functional connectivity
% ran each neuron independently
