%%%TODO'S!!!!
%{
%put all the plots in a folder of plots
%}

%% Main protocol for the experiment
%Start by running:
%test_bmi_nidaq_triggers.m

% define Animal, day and folder where to save
animal = 'NY20'; day = 'BMItest';
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
%Import: Mark Point Series, Bot half of MarkPoints

%% Update prairie view repetitions based on num neurons to stim
num_stim_neurons = max(max(holoMask));
stim_time_per_neuron = 2.1;
disp(num_stim_neurons*stim_time_per_neuron*frameRate)

%% Run HOLO STIM
%This does one neuron at a time.
clear s
HoloAcqnvsPrairie(folder, animal, day, holoMask)
%WE WANT RED+GREEN FOR THIS

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
% load the VoltageRec to check the results of holoStim
plotHoloStimTimeLock(holoActivity, voltageRec, 40, 1000) % --> To plot the result of
%ToDo: for plotting, do sliding window deltaf/f

%% Baseline acquisition
% loads the result of OnAcid
% REMOVE RED CHANNEL!!!!

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

plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)

%%
%Manually enter:
E1_base = sort([14 13 18 17], 'ascend') %JUST NEEDS GCAMP
E2_base = sort([22 3 7 5], 'ascend') %NEEDS CHROME

ensembleNeurons = [E1_base, E2_base];
plotNeuronsEnsemble(baseActivity, ensembleNeurons)

% E2_candidates = [39 45 59 37 88 6 26 46 78 48 22 20 33]
%TODO: 
%onacid from the baseline acquisition?
%if we don't do onacid on the baseline activity, we won't have any
%roi's just from gcamp

%% Calibrate Target with Baseline simulation
% select correct parameters on
% vivek_tb_test_baseline_to_calibration

% base_file = fullfile(savePath, 'BaselineOnline190513T232418.mat')
base_file = fullfile(savePath, 'BaselineOnline.mat')
exist(base_file)
n_f_file = base_file;
exist(n_f_file)
ndata = load(n_f_file);
num_base_samples = sum(~isnan(ndata.baseActivity(1,:))); 
baseline_frameRate = num_base_samples/(15*60);
A_file = fullfile(savePath, 'red.mat'); 
exist(A_file)
onacid_bool = 0

sec_per_reward_range = [80 60]; 

frames_per_reward_range = sec_per_reward_range*baseline_frameRate %[1 1.5]*60*frameRate
%multiply by frames per minute to convert
%to 

target_on_cov_bool = 0
prefix_win = 40
f0_win_bool = 1
f0_win = 2*60*ceil(frameRate)
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
    EnsembleMask = zeros(size(holoMask,1));
    for indn = 1:length(E2_base)
        auxmask = holoMask;
        auxmask(auxmask~=E2_base(indn)) = 0;
        auxmask(auxmask~=0) = indn;
        EnsembleMask = auxmask + EnsembleMask;
    end
end
figure;
imshow(EnsembleMask); 

pl = actxserver('PrairieLink.Application');
pl.Connect();
savePrairieFiles(savePath, pl, '4neurons_together')
createGplFile(savePath, EnsembleMask, posz, px, 'ensemble_')

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


%% Runs connectivity
%define where to save the file
if pl.Connected()
    pl.Disconnect();
end
savePrairieFiles(savePath, pl, 'connectivity_pre')

% create the stim train for prairie
createXmlFile(savePath, 4, 5, 0.2, 'connectivity_pre', true)
%Import: Mark Point Series, Bot half of MarkPoints
% Update prairie view repetitions based on num neurons to stim

clear s
ConnectivityAcqnvsPrairie(folder, animal, day, ensembleMask, 'PRE')

%%
%Compute vectorHolo
frameRate = 30 %baseline_frameRate
expectedLengthExperiment = 40*60*frameRate
[vectorHolo, ISI] = createVectorHolo(frameRate, expectedLengthExperiment, 20, 5, false);

%% create masks bot and image to check during experiment
createBot(savePath, x(ensembleNeurons),y(ensembleNeurons))
ensembleMask = holoMaskRedGreen;
ensembleMask(~ismember(ensembleMask,ensembleNeurons))= 0;
figure('Position', [600,300, 256, 256])
imshow(ensembleMask);

%% run Pre-training
%Change the Mark Points:
%Make 120 reps, put "Wait for Trigger" = First Reptition, Trigger
%Selection: Start with External, PFI1
clear s
baselineCalibrationFile = 'BMI_target_info_20190514T015556.mat';
vectorVTA = []
%expt_str: 
%     expt_cell = {...
%         'BMI', ...
%         'HoloVTA_pretrain', ...
%         'Holo_pretrain', ...
%         'VTA_pretrain'}; 

expt_str = 'BMI'; 
debug_bool = 0; 
debut_input = []; 
BMIAcqnvsPrairienoTrialsHoloCL_debug_enable(folder, animal, day, ...
    expt_str, baselineCalibrationFile, frameRate, vectorHolo, vectorVTA, ...
    cursor_zscore_bool, debug_bool, debug_input)
% BMIAcqnvsPrairienoTrialsHoloCL(folder, animal, day, expt_str, baselineCalibrationFile, baseline_frameRate, vectorHolo, vectorVTA, cursor_zscore_bool)

%% run BMI
%remember to set the markpoints to proper stim
BMIAcqnvsPrairienoTrials(folder, animal, day, 'BMI', baselineCalibrationFile, frameRate, vectorHolo, vectorVTA)

%% Holo stim for functional connectivity
% upload the GPL file
%% Runs connectivity
%define where to save the file
if pl.Connected()
    pl.Disconnect();
end
savePrairieFiles(savePath, pl, 'connectivity_post')

% create the stim train for prairie
createXmlFile(savePath, 4, 5, 0.2, 'connectivity_post', true)
%Import: Mark Point Series, Bot half of MarkPoints
% Update prairie view repetitions based on num neurons to stim

clear s
ConnectivityAcqnvsPrairie(folder, animal, day, ensembleMask, 'POST')

% ran each neuron independently
 