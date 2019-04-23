% Main protocol for the experiment

%% Select ChroME neurons
% Select area of interest with the 2p
% find red neurons
[holoMask, Im] = makeMasksPrairie(channel);
% it returns the mask that will be used for the holostim

% delete neurons that we don't want by position on image
holoMask = deleteMask(Im, holoMask, numArea);

%% HOLO STIM 
% creates holos with the mask of the red components as input
createGplFile(savePath, holoMask, posz)
% upload the .gpl file in the SLM
% run tseries with holoStim one neuron at a time ~2min and save the images.

%% Obtain spatial components
% obtain A from onacid, compare to red neurons and bring it to matlab
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
BaselineAcqnvsPrairie(animal, day, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

%% Baseline simulation
baseline2target(n_f_file, Acomp_file, E1_base, E2_base, frames_per_reward_range, ...
    target_on_cov_bool, prefix_win, f0_win_bool, f0_win, dff_win_bool, dff_win, save_dir)
% frames_per_reward_range must be higher than 80seconds (to keep the
% occurence of artificial vs natural higher than 80% 
% Selection of neurons


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
