%mainProt_stim_res
%I DIDNT FINISH THIS, See mainProt_grid_stim.m
%Feb 13, 2020
%measure the power out of the imaging and stim laser
%Tools > maintenance > galvo center
%(imaging laser hard shutter)
%(prairieview shutter override)


%% DEFINE PATHS

cd G:\VivekNuria\Code\HoloBMI
%DEFINE PATH_DATA: 
%
%LOAD PATHS: 
load_path = define_and_load_bmi_paths()

%SAVE PATHS: 
home_dir = 'G:\VivekNuria\Code\HoloBMI'
cd(home_dir)
env_dir = 'G:\VivekNuria\utils'

% define Animal, day and folder where to save
animal = 'NVI17'; day = 'D19';
folder = 'E:\holobmi_E\191123';
savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
redgreen_dir = fullfile(savePath, 'redgreen'); 
if ~exist(redgreen_dir, 'dir')
    mkdir(redgreen_dir);
end

path_data.load_path = load_path; 
path_data.home_dir = home_dir; %home_dir
path_data.env_dir = env_dir; %contains env files for prairie
path_data.savePath = savePath; 
path_data.im = fullfile(savePath, 'im'); %directory for imaging data
if ~exist(path_data.im, 'dir')
    mkdir(path_data.im);
end

%%
%DO: enter zoom (either 1.5 or 2)
zoom = 2; 
posz = 0;
pl = actxserver('PrairieLink.Application');
pl.Connect(); disp('Connecting to prairie...');
% Prairie variables
px = pl.PixelsPerLine();
py = pl.LinesPerFrame();
micronsPerPixel.x = str2double(pl.GetState('micronsPerPixel', 'XAxis')); 
micronsPerPixel.y = str2double(pl.GetState('micronsPerPixel', 'YAxis')); 
pl.Disconnect();
disp('Disconnected from prairie');
frameRate = 29.989;
chan_data = struct(...
    'label', 'g', ...
    'chan_idx', 2); %in RGB, G is 2nd
num_chan = length(chan_data); 
onacid_bool = false;


%%
%--------------------------------------------------------------------------
%DO:
%find FOV
%-disable the motor control!!!!
%% REMEMBER TO TURN OFF PHASE OFFSET
%% TURN OFF THE MANIPULATOR 
%% TURN OFF AUTOSCALE
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

%%
%Option 1: load images from prairie directly
option1_bool = 0; 
if option1_bool
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    disp('Connecting to prairie')
    pause(2);    
    green_im = pl.GetImage_2(2, px, py);
    red_im = pl.GetImage_2(1, px, py);
    pl.Disconnect();
    disp('Disconnected from prairie')
else
    red_path      = fullfile(redgreen_dir, 'red.tif'); 
    exist(red_path)
    green_path    = fullfile(redgreen_dir, 'green.tif'); 
%     green_path    = fullfile(redgreen_dir, 'green_mean.tif'); 
    exist(green_path)    
    green_im    = imread(green_path); 
    red_im      = imread(red_path); 
end

%%
im_summary = green_im;
im_sc_struct = struct(...
    'im', [], ...
    'minmax_perc', [], ...
    'minmax', [], ...
    'min', [], ...
    'min_perc', [], ...
    'max', [], ...
    'max_perc', []); 
num_im_sc = 0;

%Use green channel to draw ROIs
%Scale the image, in order to help see ROIs better.
%If no modification to original image needed, just run code, in 
% 'scale_im_interactive' set min_perc = 0, max_perc = 100
[im_sc_struct, num_im_sc] = scale_im_interactive(im_summary, im_sc_struct, num_im_sc);
close all;


%%
%--------------------------------------------------------------------------
%DO:
%-Input index to 'im_sc_struct'
%--------------------------------------------------------------------------
%Can input a different index to choose as the Image for choosing ROI
%Defaults to the last image in 'im_sc_struct'
im_bg = im_sc_struct(end).im; 
h = figure;
imagesc(im_bg); 
axis square
colormap('gray'); 
title('selected background image for identifying ROI'); 
%PLOT_IMAGES data:
%'plot_images' contains a set of images so user can tell if ROI selection is
%appropriate.
plot_images = struct('im', [], 'label', ''); 
plot_images(1).im = im_summary; 
plot_images(1).label = 'green mean';
plot_images(2).im = im_bg; 
plot_images(2).label = 'scaled';

%% INIT ROI_DATA
%--------------------------------------------------------------------------
%DO:
%-Confirm parameters for automatically identifying ROI
%-Run this cell
%-if template matching sucks, can set 'auto_init' to 0
%--------------------------------------------------------------------------
auto_init = 0;  %initializes roi_data using automatic cell detection: 
% Parameters for auto cell detection:
% Following were for zoom=2 on bruker soma:
% template_diam = 25; %diamter of difference of Gaussians in pixels
% thres = 0.5; %cell detection threshold as correlation coefficient
% cell_diam = 7; %CELL_DIAM is diameter used for dilation.
% finemode = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
% temmode = 0;  % 0 is for full circle (soma) 1 is for donuts (membrane)
template_diam = 15; %diamter of difference of Gaussians in pixels
thres = 0.6; %cell detection threshold as correlation coefficient
cell_diam = 14; %CELL_DIAM is diameter used for dilation.
finemode = 1; %imTemplateMatch will be used instead of normxcorr2. It will be slower.
temmode = 1;  % 0 is for full circle (soma) 1 is for donuts (membrane)
if auto_init
    %FIND ROI AUTOMATICALLY 
    [mask_intermediate, ~] = imFindCellsTM (im_bg, template_diam, thres, cell_diam, finemode, temmode);
    init_roi_mask = bwlabel(mask_intermediate);
    findCenter (init_roi_mask, im_bg);
    roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, chan_data);
else
    roi_data = init_roi_data(im_bg, num_chan, chan_data);
end

% use the time now to draw some neurons in the screen
% if result looks weird check pixel size, zoom and parameters of the
% function
% if rois are too small, remove them and draw them yourself

%%
%Visualize: 
screen_size = get(0,'ScreenSize');
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
hold on;
imagesc(roi_data.im_roi); %colormap('gray');  
axis square;
title('ROI footprint overlay in blue'); 
% scatter(roi_data.x, roi_data.y, pi*roi_data.r.^2, 'r'); 
h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
% hold on;
imagesc(roi_data.roi_mask); %colormap('gray');  
axis square;
title('ROI Mask'); 
% scatter(roi_data.x, roi_data.y, pi*roi_data.r.^2, 'r'); 
%TESTS
% roi_data = label_mask2roi_data_single_channel(im_bg, init_roi_mask, chan_data);
% [roi_ctr] = roi_bin_cell2center_radius(roi_data.roi_bin_cell);

%% Add ROI if needed
%--------------------------------------------------------------------------
%DO:
%-Run cell to manually draw additional ROI
%--------------------------------------------------------------------------
disp('Adding ROIs to image!'); 
[roi_data] = draw_roi_g_chan(plot_images, roi_data);
close all
%To Do: make it save the ROI data with each addition, so if it crashes you
%don't lose all the drawings.

%% Delete ROI if needed
%--------------------------------------------------------------------------
%DO:

%-If auto detected ROI suck, delete ROI
%--------------------------------------------------------------------------
close all;
%Delete ROI if needed: 
disp('Deleting ROIs from image!');
[roi_data] = delete_roi_2chan(plot_images, roi_data);
close all;

%% SEE ROI if needed
see_roi_data = 1; 
if see_roi_data
    %HoloMask
    holoMask = roi_data.roi_mask; 
    screen_size = get(0,'ScreenSize');
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);    
    imagesc(holoMask)
    axis square; 
    title(['holoMask num roi: ' num2str(roi_data.num_rois)]); 
    
    %im_roi:
    screen_size = get(0,'ScreenSize');
    h = figure('Position', [screen_size(3)/2 1 screen_size(3)/2 screen_size(4)]);
    % hold on;
    imagesc(roi_data.im_roi); %colormap('gray');  
    axis square;
    title(['ROI footprint overlay in blue.  Num ROI: ' num2str(roi_data.num_rois)]);
end

%% Save roi_data
roi_data_file = fullfile(savePath, 'roi_data.mat'); 
roi_mask = roi_data.roi_mask;
save(roi_data_file, 'roi_mask', 'plot_images', 'im_sc_struct', 'roi_data'); 

%% prepare SEQUENTIAL HOLO STIM of individual neurons
%Load environment to prairie
% load environment
env_path = fullfile(env_dir, "Tseries_VivekNuria_holo_all.env")
exist(env_path)
pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + env_path;
pl.SendScriptCommands(loadCommand);
pl.Disconnect()

%% %Create gpl, xml files for individual points
%DEFINE GPL FILE
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', 100, ...
    'SpiralSize', 0.7, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals

markpoints_data = repmat(init_markpoints, [roi_data.num_rois 1]); 

%INPUT: USE A FIXED SPIRAL SIZE?
use_fixed_size_bool = 0;
fixed_spiral_size_um = 15
spiral_size_gain = 1 %set this to scale the size of each mark points circle
for roi_i = 1:roi_data.num_rois
    if  used_fixed_size_bool
        spiral_size_um = fixed_spiral_size_um
    else
        spiral_size_um  = 2*roi_data.r(roi_i)*micronsPerPixel.x*spiral_size_gain;
    end
    spiral_size     = spiral_size_conversion*spiral_size_um;
    markpoints_data(roi_i).SpiralSize = spiral_size; %double(ceil(10*spiral_size)/10);
end
% creates holos
createGplFile_v2(savePath, markpoints_data, roi_data.x, roi_data.y, posz, roi_data.r, px, zoom)

%% XML: Sequential Single Cell Stim
xml_seq_path = fullfile(savePath, 'seq_single_stim.xml'); 
%{
Params Summary:
-num_sequences

-Initial Delay: time (ms) between each stimulation
-power
-numSpirals
-Repetitions
-Iter - number of iterations
-IterDelay - time beteen iterations
-InterPointDelay

%}
numberNeurons= roi_data.num_rois; %max(max(holoMask));
% numberNeurons=max(max(holoMask));
num_sequences = 1;
stim_sequence = repmat(1:numberNeurons, [1 num_sequences]); 

power_conversion = 0.004; %0.2 -> 50, 0.4->100

seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
%
InitialDelay = 2000; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,numberNeurons);
%
power = 40;
power_converted = power*power_conversion;
seq_stim_params.PowerVector = power_converted*ones(1,numberNeurons);
%2
Duration = 30;
seq_stim_params.DurationVector = Duration*ones(1,numberNeurons);
% 
numSpirals = 10;
seq_stim_params.SpiralVector = numSpirals*ones(1,numberNeurons);
%
Repetitions = 1;
seq_stim_params.RepetitionsVector = Repetitions*ones(1,numberNeurons);
%Darcy sometimes recommends increasing 'Repetitions' and decreasing
%Duration.  This changes the distribution of the spirals in time over the
%cell.
%
InterPointDelay =  0.12;
seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,numberNeurons); 

% seq_stim_params
createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, stim_sequence);

% Update prairie view repetitions based on num neurons to stim
stim_time_per_neuron = InitialDelay/1000+InterPointDelay;
num_reps_seq_stim = ceil(numberNeurons*stim_time_per_neuron*frameRate);
len_seq_stim = numberNeurons*stim_time_per_neuron/60;

disp(['Number of Repetitions in PrairieView: ' num2str(num_reps_seq_stim)])
disp(['Stim time per neuron (s): ' num2str(stim_time_per_neuron)]); 
disp(['Num neurons: ' num2str(numberNeurons)]); 
disp(['Length (min): ' num2str(len_seq_stim)])

% if position of stim cells looks different, "smaller/bigger" check the
% pixel size
%--------------------------------------------------------------------------
%DO: 
% upload .gpl in MarkPoints (Top half)
% upload .xml in MarkPoints (Bot half)
% update T-series repetitions in Prairie View with above number
% Make sure Voltage Recording has all channels enabled
% Make sure you turn on the laser power and pmt's
% while running paint the neurons
% TODO: automate uploading
%--------------------------------------------------------------------------

%% Run HOLO STIM to check stim-able neurons
%This stims one neuron at a time.
%--------------------------------------------------------------------------
%DO: 
%1) Do live scan check
%--------------------------------------------------------------------------

clear s
expt_str = 'holostim_seq'; %previously 'holostim' 
mask = roi_data.roi_mask;
expectedLengthExperiment = ceil(num_reps_seq_stim*1.5); 
HoloAcqnvsPrairie_v2(path_data, expt_str, mask, expectedLengthExperiment)
% HoloAcqnvsPrairie(folder, animal, day, holoMask)
%TODO: make this closed loop, and wait for the neurons to be inactive
%before stimming them.
%{
%--------------------------------------------------------------------------
%DO: 
%Image-block ripping utility: Convert the holostim acqn (2 files)
%Load holostim_seqX.mat
%Load voltage recording for plotting.  
%   Use 'Import Data' in matlab
%   Import csv to matlab: output type is Numeric Matrix.  
%   Name it "voltageRec")
%--------------------------------------------------------------------------
%}

%%
%(Image-Block Ripping Utility) Convert holostim file with bruker converter 
% load the VoltageRec to check the results of holoStim
min_duration = 40; %stims can't occur within this number of samples of voltageRec
plot_win = 1000; 
plotHoloStimTimeLock(holoActivity, voltageRec, min_duration, plot_win)
%ToDo: allow us to select the idxs of neurons to plot
%ToDo: for plotting, do sliding window deltaf/
%ToDo: a version that just plots each individual neuron, we type 'Y' or 'N'
%to make it a candidate

%%
%Choose a neuron to stim: 
E2_ens = [1];

%% Holo stim of Ensemble neurons
% Make GPL (points), BOT (measure activity)
% -select markpoints_data fo
close all
sel_idxs = unique(E2_ens); 
[sel_roi_data, sel_idxs] = select_roi_data(roi_data, sel_idxs); 
sel_markpoints_data = markpoints_data(sel_idxs); 

%--
%GPL for Stim Ensemble
% savePrairieFiles(savePath, pl, 'GPL_candidates_')
% creates holos
gpl_candidates_path = createGplFile_v2(savePath, sel_markpoints_data, sel_roi_data.x, sel_roi_data.y, posz, sel_roi_data.r, px, zoom, 'GPL_ens_');
% createGplFile(savePath, StimMask, posz, px, 'ensemble_')
%BOT
bot_candidates_path = fullfile(savePath, 'BOT_ens.cfg'); 
createBot_v2(bot_candidates_path, sel_roi_data.x, sel_roi_data.y, sel_roi_data.r)
% createBot(savePath, x(E2_base),y(E2_base))
%NOTE: 
%If it can't be loaded in prairie because prairie says the file is in use
%by another program, it can be because matlab didn't release the file.
%Re run, saving to a different filenaem

%XML:
numberNeurons = length(E2_ens);
ens_default_params.numberNeurons = numberNeurons;
%SPECIFY:
ens_default_params.PowerVector          = 30*...
    power_conversion*ones(1,numberNeurons);
ens_default_params.DurationVector       = 20*...
    ones(1,numberNeurons);
ens_default_params.RepetitionsVector    = 1*...
    ones(1,numberNeurons);
ens_default_params.SpiralVector         = 10*...
    ones(1,numberNeurons);
%Don't Change:
ens_default_params.InitialDelayVector   = 0*...
    ones(1,numberNeurons);
ens_default_params.UncagingLaser        = "Monaco"; 
ens_default_params.AllPointsAtOnce      = "False"
ens_default_params.Iter                 = 1; %how many times to go through and stim each cell.
ens_default_params.IterDelay            = 0; %Time (ms) between iterations
InterPointDelay                         =  0.12;
ens_default_params.InterPointDelayVector = InterPointDelay*ones(1,numberNeurons);
xml_candidates_path = fullfile(savePath, 'XML_ens.xml'); 
createXmlFile_basic(xml_candidates_path, ens_default_params);

pl = actxserver('PrairieLink.Application');
pl.Connect();
loadCommand = "-tsl " + fullfile(path_data.env_dir, "Tseries_VivekNuria_holo_4.env");
pl.SendScriptCommands(loadCommand);
pl.Disconnect();

% Check neuron size, sometimes they are genourmous sometimes they are tiny.
% Improvise, adapt, overcome

%--------------------------------------------------------------------------
%DO: 

%1) upload the GPL file
%2) upload XML file
%3) Click BOT in Image Window, load BOT.cfg, only display ROI of interest
%4) Run BOT, and adjust (duration, power, repetitions) of each neuron's stim
%5) Note down below the params for each neuron

%%
% (sel_roi_data) target x,y,r + {angles, radii} to grid of points
% grid of points: vector of x, y
% row,col mapping to x,y? 

%%
%{
%DO:
%In the GUI Mark Points Series, remove bad neurons 
%   1) from Points/Groups, 
    2) from Mark Point Series
    3) from BOT 
%For chosen neurons, Enter Stim Parameters you like in GUI.
EXPORT each of the above, with the following names: 
gpl_path = fullfile(savePath, 'GPL_ens.gpl');
xml_path = fullfile(savePath, 'XML_ens.xml');
bot_path = fullfile(savePath, 'BOT_ens.cfg');
%}

%% --------------------------------------------------------------------------
if ~onacid_bool
    AComp = 0;
else
    load(roi_data_file);
end
% Baseline environment already removes MARKPOINTS and set the reps to 27000
holoMask = roi_data.roi_mask;
[base_mat_path, base_dat_path] = ...
    BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMask, task_settings);
% BaselineAcqnvsPrairie(folder, animal, day, AComp, holoMaskRedGreen, onacid_bool, frameRate);
% saves in [savePath, 'baselineActivity.dat'] the activity of all the
% neurons of the mask (Acomp+red)
% saves in baseOnline.mat the baseline activityim_bg

%--------------------------------------------------------------------------
%D0:
%0) Abort T-series (cuz of voltage recording)
%1) pyctrl stop
%2) load cells stop
%3) video stop
%4) Drag load cell data to folder
%5) Drag video to folder
%--------------------------------------------------------------------------

%% Plot the baseline activity: 
% plots neurons so we can select which ones we like the most 

base_file = base_mat_path; 
%Can manually enter a previous path: 
% base_file = fullfile(savePath, 'BaselineOnline190526T113422.mat')
if onacid_bool
    totalneurons = min(size(AComp,2), 20);
else
    totalneurons = max(max(roi_data.num_rois));
    CComp = [];
    YrA = []; 
end
load(base_file); 
num_neurons_plot = min(30, totalneurons); %totalneurons
% plotNeuronsBaseline(baseActivity, CComp, YrA, totalneurons)
plotNeuronsBaseline(baseActivity, CComp, YrA, num_neurons_plot)

%%
