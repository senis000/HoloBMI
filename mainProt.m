% Main protocol for the experiment

% Select area of interest with the 2p

% find red neurons
[holoMask, Im] = makeMasksPrairie(channel);
% it returns the mask that will be used for the holostim

% delete neurons that we don't want by position on image
holoMask = deleteMask(Im, holoMask, numArea);

% HOLO STIM with the mask of the red components as input!!!!
% create the .gpl file with the redmask
% update the .gpl file in the SLM


% run the holostim ~2min and save the images.
% TODO:


% obtain A from onacid and bring it to matlab
% "obtain_componenents.py" -> 
%         'redLabel' : label as True for components marked as red,
%         'indRed' : index of such neurons,
%         'AComp' : sparse matrix with the spatial filters,
%         'com' : position of the neurons
% it also saves figures for sanity check

% runs the baseline acquisiton
BaselineAcqnvsPrairie(animal, day, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activity

% Baseline simulation
% Vivek
% obtains the value of T1 that will act as threshold

% run BMI
% 2 versions 
BMIAcqnvsPrairie(animal, day, E1, E2, T1, frameRate)  %not updated
BMIAcqnvsPrairienoTrials(animal, day, E1, E2, T1, frameRate)
% should be the same function with or without trials

% VTA Stim during BMI
% TODO check port, pin and duration of trigger

% HOLO Stim during BMI
% TODO
