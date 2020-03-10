
%%
%--------------------------------------------------------------------------
%BEFORE ANIMAL IN BOX:
%DO:
% Hook up BNCs: 
% 1) BMI solenoid, AI5
% 2) Monaco Trig, AI6
% 3) Frame Trig, AI7
% 4) Holo Trig PFI1
% 5) run: %main_test_190923

%%
% %Load example roi_data information for reference: 
% roi_data_path = fullfile('Z:\Vivek\holoBMI_2round\191119\NVI20\D15', 'roi_data.mat'); 
% d = load(roi_data_path)
% workspace_path = fullfile('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-13', 'workspace.mat'); 
% w = load(workspace_path)

%%
% C:\ProgramData\Bruker Fluorescence Microscopy\Prairie View\5.4.64.700\Logs

%% DEFINE PATHS

cd G:\VivekNuria\Code\HoloBMI
%DEFINE PATH_DATA: 

% define Animal, day and folder where to save
%Animals: NVI20, NVI21
animal = 'NVI20'; day = 'D0';
folder = 'E:\holobmi_E\200309p2' %'E:\holobmi_E\200305';
savePath = fullfile(folder, animal,  day);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
redgreen_dir = fullfile(savePath, 'redgreen'); 
if ~exist(redgreen_dir, 'dir')
    mkdir(redgreen_dir);
end

%LOAD PATHS: 
load_path = define_and_load_bmi_paths()
%SAVE PATHS: 
home_dir = 'G:\VivekNuria\Code\HoloBMI'
cd(home_dir)
env_dir = 'G:\VivekNuria\utils'

path_data.load_path = load_path; 
path_data.home_dir = home_dir; %home_dir
path_data.env_dir = env_dir; %contains env files for prairie
path_data.savePath = savePath; 
path_data.im = fullfile(savePath, 'im'); %directory for imaging data
if ~exist(path_data.im, 'dir')
    mkdir(path_data.im);
end

%TODO streamlin
addpath('G:\VivekNuria\Code\HoloBMI\grid_stim')

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
%find FOV (put the target cell in the middle)
%
%-disable the motor control!!!!
% REMEMBER TO TURN OFF PHASE OFFSET
% TURN OFF THE MANIPULATOR 
% TURN OFF AUTOSCALE
%-ACTIVATE RED+GREEN channels
%-collect 4000 frame video, save it in im with name "redgreen"
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
%     red_path      = fullfile(redgreen_dir, 'red.tif'); 
%     exist(red_path)
    green_path    = fullfile(redgreen_dir, 'green.tif'); 
%     green_path    = fullfile(redgreen_dir, 'green_mean.tif'); 
    exist(green_path)    
    green_im    = imread(green_path); 
%     red_im      = imread(red_path); 
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

%%
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

%%
%Choose the target roi: 
% First show the ROI: 
h = figure; imagesc(roi_data.roi_mask); axis square; colorbar

%% Choose the roi: 
target_idx = 1; 
[target_roi_data, target_idx] = select_roi_data(roi_data, target_idx)

%%

%%
%--------------------------------------------------------------------------
%STIM RELATED CODE HERE
%--------------------------------------------------------------------------
%Define target location of stim
%Default to center of image: 
%Can manually select x, y based on roi_data
target.x = target_roi_data.x %dim/2
target.y = target_roi_data.y %dim/2
target.r = 10  %not super necessary for stim grid points.  needed for mark points generation

%--------------------------------------------------------------------------
%BOT:
bot_target_path = fullfile(savePath, 'BOT_target.cfg'); 
createBot_v2(bot_target_path, target.x, target.y, target.r)

%--------------------------------------------------------------------------
%GPL:
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
initSpiralSize_um   = 14; 
initSpiralSize = spiral_size_conversion*initSpiralSize_um;
stim_duration       = 30;  %(ms)
init_markpoints     = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', stim_duration, ...
    'SpiralSize', initSpiralSize, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals
markpoints_data = repmat(init_markpoints, [1 1]); 
createGplFile_v2(savePath, markpoints_data, target.x, target.y, posz, NaN, px, zoom, 'target_');

%--------------------------------------------------------------------------
%XML: 
num_stims    = 1;
stim_sequence = ones(1,num_stims); 

xml_seq_path        = fullfile(savePath, 'target_stim.xml'); 
time_between_stims  = 0; %(ms)
power_conversion = 0.004; %0.2 -> 50, 0.4->100
seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
InitialDelay = time_between_stims; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,num_stims);
%
power = 30;
power_converted = power*power_conversion;
seq_stim_params.PowerVector = power_converted*ones(1,num_stims);
%2
Duration = stim_duration;
seq_stim_params.DurationVector = Duration*ones(1,num_stims);
% 
numSpirals = 10;
seq_stim_params.SpiralVector = numSpirals*ones(1,num_stims);
%
Repetitions = 1;
seq_stim_params.RepetitionsVector = Repetitions*ones(1,num_stims);
%Darcy sometimes recommends increasing 'Repetitions' and decreasing
%Duration.  This changes the distribution of the spirals in time over the
%cell.
%
InterPointDelay =  0.12;
seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,num_stims); 
% seq_stim_params
%--------------------------------------------------------------------------
createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, stim_sequence);

%%
%Stim each ROI two times, with 5 sec separation: 
%--------------------------------------------------------------------------
%Create BOT: 
bot_roi_path = fullfile(savePath, 'BOT_roi.cfg'); 
createBot_v2(bot_roi_path, roi_data.x, roi_data.y, roi_data.r)
num_roi = length(roi_data.x); 

%%
%--------------------------------------------------------------------------
%GPL:
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
initSpiralSize_um   = 14; 
initSpiralSize = spiral_size_conversion*initSpiralSize_um;
stim_duration       = 30;  %(ms)
init_markpoints     = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', stim_duration, ...
    'SpiralSize', initSpiralSize, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals
markpoints_data = repmat(init_markpoints, [num_roi 1]); 
roi_gpl_prefix = 'roi_'; 
createGplFile_v2(savePath, markpoints_data, roi_data.x, roi_data.y, posz, roi_data.r, px, zoom, roi_gpl_prefix);

%%
%--------------------------------------------------------------------------
%XML: 
stim_per_roi = 3; 
numberNeurons           = roi_data.num_rois
roi_stim_sequence       = repmat(1:roi_data.num_rois, [1 stim_per_roi]); %repeat a stim of each ROI three times
num_stims               = length(roi_stim_sequence);

xml_seq_path        = fullfile(savePath, 'roi_stim.xml'); 
time_between_stims  = 5000; %(ms)
power_conversion = 0.004; %0.2 -> 50, 0.4->100
seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
InitialDelay = time_between_stims; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,num_stims);
%
power = 30;
power_converted = power*power_conversion;
seq_stim_params.PowerVector = power_converted*ones(1,num_stims);
%2
Duration = stim_duration;
seq_stim_params.DurationVector = Duration*ones(1,num_stims);
% 
numSpirals = 10;
seq_stim_params.SpiralVector = numSpirals*ones(1,num_stims);
%
Repetitions = 1;
seq_stim_params.RepetitionsVector = Repetitions*ones(1,num_stims);
%Darcy sometimes recommends increasing 'Repetitions' and decreasing
%Duration.  This changes the distribution of the spirals in time over the
%cell.
%
InterPointDelay =  0.12;
seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,num_stims); 
% seq_stim_params
%--------------------------------------------------------------------------
xml_seq_path
createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, roi_stim_sequence);

% Update prairie view repetitions based on num neurons to stim
stim_time_per_neuron = InitialDelay/1000+InterPointDelay;
num_reps_seq_stim = ceil(numberNeurons*stim_time_per_neuron*frameRate*stim_per_roi);
len_seq_stim = numberNeurons*stim_time_per_neuron*stim_per_roi/60;

disp(['Number of Repetitions in PrairieView: ' num2str(num_reps_seq_stim)])
disp(['Stim time per neuron (s): ' num2str(stim_time_per_neuron)]); 
disp(['Num neurons: ' num2str(numberNeurons)]); 
disp(['Length (min): ' num2str(len_seq_stim)])


%--------------------------------------------------------------------------
%% Run HOLO STIM to check stim-able neurons
%This stims one neuron at a time.
%--------------------------------------------------------------------------
%DO: 
%--------------------------------------------------------------------------
%DO: 
% upload .gpl in MarkPoints (Top half)
% upload .xml in MarkPoints (Bot half)
% update T-series repetitions in Prairie View with above number
% Make sure Voltage Recording has all channels enabled
% Make sure you turn on the laser power and pmt's
% while running paint the neurons
% TODO: automate uploading
%1) Do live scan check
%--------------------------------------------------------------------------

clear s
expt_str = 'roi_stim_seq'; %previously 'holostim' 
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
% load the VoltageRec to check the results of holoStim (voltageRec)
%TO DO: FIX THIS: 
% min_duration = 40; %stims can't occur within this number of samples of voltageRec
% plot_win = 1000; 
% plotHoloStimTimeLock(holoActivity, voltageRec, min_duration, plot_win)
%
%ToDo: allow us to select the idxs of neurons to plot
%ToDo: for plotting, do sliding window deltaf/
%ToDo: a version that just plots each individual neuron, we type 'Y' or 'N'
%
%to make it a candidate


%%
%--------------------------------------------------------------------------
%Coarse Grid: 
%1) define grid in pixel space
% micronsPerPixel.x = .7266
% micronsPerPixel.y = .7266
%USER INPUT: 
dim = 512
fine_grid = 0

if fine_grid
    x_extent_microns = [-20 20]; 
    x_extent = round(x_extent_microns/micronsPerPixel.x)
    y_extent_microns = [-20 20]; 
    y_extent = round(y_extent_microns/micronsPerPixel.y)
    
    step_microns = 5
    x_step = round(step_microns/micronsPerPixel.x);
    y_step = round(step_microns/micronsPerPixel.y);    
    
else
    
    x_extent_microns = [-60 60]; 
    x_extent = round(x_extent_microns/micronsPerPixel.x)
    y_extent_microns = [-60 60]; 
    y_extent = round(y_extent_microns/micronsPerPixel.y)
    
    step_microns = 10
    x_step = round(step_microns/micronsPerPixel.x);
    y_step = round(step_microns/micronsPerPixel.y);    
    
end
[x_mesh_flat, y_mesh_flat, ctr_idx, num_grid_pts, x_mesh, y_mesh] = ...
    gen_square_grid(target, dim, x_extent, y_extent, x_step, y_step);
num_grid_pts

% [x_mesh_flat, y_mesh_flat, ctr_idx, num_grid_pts] = ...
%     gen_L_grid(target, dim, x_extent, y_extent, x_step, y_step);

% %Directly set extent in pixel space: 
% x_extent = [-dim dim]/4
% y_extent = [-dim dim]/4

h = figure; imagesc(im_bg); colormap('gray'); axis square; 
hold on
scatter(x_mesh_flat, y_mesh_flat, 'x'); 
scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 'r', 'x'); 
xlim([1 dim])
ylim([1 dim])
axis square
%Reminder, mapping from x,y to image space: 
%x = column, y = row

%
%--------------------------------------------------------------------------
%2) map grid to stim points (GPL file): 
%DEFINE GPL FILE
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
initSpiralSize_um = 14; 
initSpiralSize = spiral_size_conversion*initSpiralSize_um;

stim_duration       = 30;  %(ms)

init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', stim_duration, ...
    'SpiralSize', initSpiralSize, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals

markpoints_data = repmat(init_markpoints, [num_grid_pts 1]); 
grid_gpl_prefix = 'coarse_grid_'; 
createGplFile_v2(savePath, markpoints_data, x_mesh_flat, y_mesh_flat, posz, NaN, px, zoom, grid_gpl_prefix);

% If we want each grid point to be a different radius stim: 
% use_fixed_size_bool = 0;
% fixed_spiral_size_um = 15
% spiral_size_gain = 1 %set this to scale the size of each mark points circle
% for roi_i = 1:roi_data.num_rois
%     if  used_fixed_size_bool
%         spiral_size_um = fixed_spiral_size_um
%     else
%         spiral_size_um  = 2*roi_data.r(roi_i)*micronsPerPixel.x*spiral_size_gain;
%     end
%     spiral_size     = spiral_size_conversion*spiral_size_um;
%     markpoints_data(roi_i).SpiralSize = spiral_size; %double(ceil(10*spiral_size)/10);
% end

%%
%--------------------------------------------------------------------------
%3) Make a sequence of stimulations: 
%Every 'num_grid_stim_per_target_stim', stim the target point.   
num_reps = 2; %how many times should we stim each grid point
num_grid_stim_per_target_stim = 20;

num_stims_expected = num_reps*(num_grid_pts + ceil(num_grid_pts/num_grid_stim_per_target_stim))

rand_order = repmat(randperm(num_grid_pts), [1 num_reps]); 
stim_sequence = [ctr_idx ctr_idx ctr_idx]

for rep_i = 1:num_reps
    for i = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
        sel_from_rand_order = ...
            (1:num_grid_stim_per_target_stim)+(i-1)*num_grid_stim_per_target_stim;   
        sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
        grid_idxs = rand_order(sel_from_rand_order); 
        for j = 1:length(grid_idxs)
            stim_sequence = [stim_sequence grid_idxs(j)];
        end
        stim_sequence = [stim_sequence ctr_idx];
    end
end
% h = figure;
% plot(stim_sequence)
%
%Old way of generating the sequence: 
% for i = 1:num_reps
%     rand_order = randperm(num_grid_pts);
%     for j = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
%         sel_from_rand_order = ...
%             (1:num_grid_stim_per_target_stim)+(j-1)*num_grid_stim_per_target_stim;
%         %Remove invalid idxs
%         sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
%         grid_idxs = rand_order(sel_from_rand_order); 
%         stim_sequence = [stim_sequence ctr_idx grid_idxs];
%     end
% end
num_stims = length(stim_sequence)
expected_duration = num_stims*(time_between_stims/1000)/60
%Confirmation:
h= figure;
hist(stim_sequence, 1:num_grid_pts)
title(['num stims per grid point, ctr idx: ' num2str(ctr_idx)])
disp('center point:')
ctr_idx

% XML: Sequential Single Cell Stim

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
%INPUT:
stim_duration       = 50; 
time_between_stims  = 8000; %(ms)
expectedDuration    = (stim_duration + time_between_stims)*num_stims/(60*1000) %min 

%
xml_seq_path = fullfile(savePath, 'coarse_grid_stim.xml'); 
power_conversion = 0.004; %0.2 -> 50, 0.4->100

seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
%
InitialDelay = time_between_stims; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,num_stims);
%
power = 30;
power_converted = power*power_conversion;
seq_stim_params.PowerVector = power_converted*ones(1,num_stims);
%2
Duration = stim_duration;
seq_stim_params.DurationVector = Duration*ones(1,num_stims);
% 
numSpirals = 10;
seq_stim_params.SpiralVector = numSpirals*ones(1,num_stims);
%
Repetitions = 1;
seq_stim_params.RepetitionsVector = Repetitions*ones(1,num_stims);
%Darcy sometimes recommends increasing 'Repetitions' and decreasing
%Duration.  This changes the distribution of the spirals in time over the
%cell.
%
InterPointDelay =  0.12;
seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,num_stims); 

% seq_stim_params
%--------------------------------------------------------------------------
createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, stim_sequence);
%--------------------------------------------------------------------------
% Update prairie view repetitions based on num neurons to stim
stim_time_per_neuron = InitialDelay/1000+InterPointDelay;
len_seq_stim = numberNeurons*stim_time_per_neuron/60;

num_reps_seq_stim = ceil(expectedDuration*60*frameRate);

disp(['Number of Repetitions in PrairieView: ' num2str(num_reps_seq_stim)])
disp('Update PrairieView to have at least this many reps'); 
disp(['Length (min): ' num2str(expectedDuration)])

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
%Save the matlab workspace data!

%%
% clear s
% expt_str = 'gridstim'; %previously 'holostim' 
% % expt_str = 'test_gridstim'; %previously 'holostim' 
% mask = roi_data.roi_mask;
% expectedLengthExperiment = ceil(num_reps_seq_stim*1.1); 
% HoloAcqnvsPrairie_v2(path_data, expt_str, mask, expectedLengthExperiment)









%%
%--------------------------------------------------------------------------
%Fine Grid: 
%1) define grid in pixel space
% micronsPerPixel.x = .7266
% micronsPerPixel.y = .7266
%USER INPUT: 
dim = 512
fine_grid = 1

if fine_grid
    x_extent_microns = [-20 20]; 
    x_extent = round(x_extent_microns/micronsPerPixel.x)
    y_extent_microns = [-20 20]; 
    y_extent = round(y_extent_microns/micronsPerPixel.y)
    
    step_microns = 5
    x_step = round(step_microns/micronsPerPixel.x);
    y_step = round(step_microns/micronsPerPixel.y);    
    
else
    
    x_extent_microns = [-60 60]; 
    x_extent = round(x_extent_microns/micronsPerPixel.x)
    y_extent_microns = [-60 60]; 
    y_extent = round(y_extent_microns/micronsPerPixel.y)
    
    step_microns = 10
    x_step = round(step_microns/micronsPerPixel.x);
    y_step = round(step_microns/micronsPerPixel.y);    
    
end
[x_mesh_flat, y_mesh_flat, ctr_idx, num_grid_pts, x_mesh, y_mesh] = ...
    gen_square_grid(target, dim, x_extent, y_extent, x_step, y_step);
num_grid_pts

% [x_mesh_flat, y_mesh_flat, ctr_idx, num_grid_pts] = ...
%     gen_L_grid(target, dim, x_extent, y_extent, x_step, y_step);

% %Directly set extent in pixel space: 
% x_extent = [-dim dim]/4
% y_extent = [-dim dim]/4

h = figure; imagesc(im_bg); colormap('gray'); axis square; 
hold on
scatter(x_mesh_flat, y_mesh_flat, 'x'); 
scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 'r', 'x'); 
xlim([1 dim])
ylim([1 dim])
axis square
%Reminder, mapping from x,y to image space: 
%x = column, y = row

%
%--------------------------------------------------------------------------
%2) map grid to stim points (GPL file): 
%DEFINE GPL FILE
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
initSpiralSize_um = 14; 
initSpiralSize = spiral_size_conversion*initSpiralSize_um;

stim_duration       = 30;  %(ms)

init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', stim_duration, ...
    'SpiralSize', initSpiralSize, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals

markpoints_data = repmat(init_markpoints, [num_grid_pts 1]); 
grid_gpl_prefix = 'fine_grid_'; 
createGplFile_v2(savePath, markpoints_data, x_mesh_flat, y_mesh_flat, posz, NaN, px, zoom, grid_gpl_prefix);

% If we want each grid point to be a different radius stim: 
% use_fixed_size_bool = 0;
% fixed_spiral_size_um = 15
% spiral_size_gain = 1 %set this to scale the size of each mark points circle
% for roi_i = 1:roi_data.num_rois
%     if  used_fixed_size_bool
%         spiral_size_um = fixed_spiral_size_um
%     else
%         spiral_size_um  = 2*roi_data.r(roi_i)*micronsPerPixel.x*spiral_size_gain;
%     end
%     spiral_size     = spiral_size_conversion*spiral_size_um;
%     markpoints_data(roi_i).SpiralSize = spiral_size; %double(ceil(10*spiral_size)/10);
% end

%%
%--------------------------------------------------------------------------
%3) Make a sequence of stimulations: 
%Every 'num_grid_stim_per_target_stim', stim the target point.   
num_reps = 2; %how many times should we stim each grid point
num_grid_stim_per_target_stim = 20;

num_stims_expected = num_reps*(num_grid_pts + ceil(num_grid_pts/num_grid_stim_per_target_stim))

rand_order = repmat(randperm(num_grid_pts), [1 num_reps]); 
stim_sequence = [ctr_idx ctr_idx ctr_idx]

for rep_i = 1:num_reps
    for i = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
        sel_from_rand_order = ...
            (1:num_grid_stim_per_target_stim)+(i-1)*num_grid_stim_per_target_stim;   
        sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
        grid_idxs = rand_order(sel_from_rand_order); 
        for j = 1:length(grid_idxs)
            stim_sequence = [stim_sequence grid_idxs(j)];
        end
        stim_sequence = [stim_sequence ctr_idx];
    end
end
% h = figure;
% plot(stim_sequence)
%
%Old way of generating the sequence: 
% for i = 1:num_reps
%     rand_order = randperm(num_grid_pts);
%     for j = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
%         sel_from_rand_order = ...
%             (1:num_grid_stim_per_target_stim)+(j-1)*num_grid_stim_per_target_stim;
%         %Remove invalid idxs
%         sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
%         grid_idxs = rand_order(sel_from_rand_order); 
%         stim_sequence = [stim_sequence ctr_idx grid_idxs];
%     end
% end
num_stims = length(stim_sequence)
expected_duration = num_stims*(time_between_stims/1000)/60
%Confirmation:
h= figure;
hist(stim_sequence, 1:num_grid_pts)
title(['num stims per grid point, ctr idx: ' num2str(ctr_idx)])
disp('center point:')
ctr_idx

% XML: Sequential Single Cell Stim

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
%INPUT:
stim_duration       = 50; 
time_between_stims  = 8000; %(ms)
expectedDuration    = (stim_duration + time_between_stims)*num_stims/(60*1000) %min 

%
xml_seq_path = fullfile(savePath, 'fine_grid_stim.xml'); 
power_conversion = 0.004; %0.2 -> 50, 0.4->100

seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
%
InitialDelay = time_between_stims; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,num_stims);
%
power = 30;
power_converted = power*power_conversion;
seq_stim_params.PowerVector = power_converted*ones(1,num_stims);
%2
Duration = stim_duration;
seq_stim_params.DurationVector = Duration*ones(1,num_stims);
% 
numSpirals = 10;
seq_stim_params.SpiralVector = numSpirals*ones(1,num_stims);
%
Repetitions = 1;
seq_stim_params.RepetitionsVector = Repetitions*ones(1,num_stims);
%Darcy sometimes recommends increasing 'Repetitions' and decreasing
%Duration.  This changes the distribution of the spirals in time over the
%cell.
%
InterPointDelay =  0.12;
seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,num_stims); 

% seq_stim_params
%--------------------------------------------------------------------------
createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, stim_sequence);
%--------------------------------------------------------------------------
% Update prairie view repetitions based on num neurons to stim
stim_time_per_neuron = InitialDelay/1000+InterPointDelay;
len_seq_stim = numberNeurons*stim_time_per_neuron/60;

num_reps_seq_stim = ceil(expectedDuration*60*frameRate);

disp(['Number of Repetitions in PrairieView: ' num2str(num_reps_seq_stim)])
disp('Update PrairieView to have at least this many reps'); 
disp(['Length (min): ' num2str(expectedDuration)])

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

%%








% %%
% %--------------------------------------------------------------------------
% %3) Make a sequence of stimulations: 
% %Every 'num_grid_stim_per_target_stim', stim the target point. 
% num_select = 30; 
% num_reps = 2; %how many times should we stim each grid point
% num_grid_stim_per_target_stim = 20;
% 
% num_stims_expected = num_reps*(num_grid_pts + ceil(num_grid_pts/num_grid_stim_per_target_stim))
% 
% rand_order = repmat(randperm(num_grid_pts), [1 num_reps]); 
% stim_sequence = [ctr_idx]
% 
% for rep_i = 1:num_reps
%     for i = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
%         sel_from_rand_order = ...
%             (1:num_grid_stim_per_target_stim)+(i-1)*num_grid_stim_per_target_stim;   
%         sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
%         grid_idxs = rand_order(sel_from_rand_order); 
%         for j = 1:length(grid_idxs)
%             stim_sequence = [stim_sequence grid_idxs(j)];
%         end
%         stim_sequence = [stim_sequence ctr_idx];
%     end
% end
% stim_sequence = [ctr_idx ctr_idx ctr_idx stim_sequence]
% stim_sequence = stim_sequence(1:num_select); 
% % h = figure;
% % plot(stim_sequence)
% %
% %Old way of generating the sequence: 
% % for i = 1:num_reps
% %     rand_order = randperm(num_grid_pts);
% %     for j = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
% %         sel_from_rand_order = ...
% %             (1:num_grid_stim_per_target_stim)+(j-1)*num_grid_stim_per_target_stim;
% %         %Remove invalid idxs
% %         sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
% %         grid_idxs = rand_order(sel_from_rand_order); 
% %         stim_sequence = [stim_sequence ctr_idx grid_idxs];
% %     end
% % end
% num_stims = length(stim_sequence)
% expected_duration = num_stims*(time_between_stims/1000)/60
% %Confirmation:
% h= figure;
% hist(stim_sequence, 1:num_grid_pts)
% title(['num stims per grid point, ctr idx: ' num2str(ctr_idx)])
% disp('center point:')
% ctr_idx
% 
% % XML: Sequential Single Cell Stim
% 
% %{
% Params Summary:
% -num_sequences
% 
% -Initial Delay: time (ms) between each stimulation
% -power
% -numSpirals
% -Repetitions
% -Iter - number of iterations
% -IterDelay - time beteen iterations
% -InterPointDelay
% 
% %}
% %INPUT:
% stim_duration       = 30; 
% time_between_stims  = 5000; %(ms)
% expectedDuration    = (stim_duration + time_between_stims)*num_stims/(60*1000) %min 
% 
% %
% xml_seq_path = fullfile(savePath, 'grid_stim_test_v3.xml'); 
% power_conversion = 0.004; %0.2 -> 50, 0.4->100
% 
% seq_stim_params.UncagingLaser = "Monaco"; 
% seq_stim_params.AllPointsAtOnce = "False"
% seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
% seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
% %
% InitialDelay = time_between_stims; %(ms) time bw stim delivery
% seq_stim_params.InitialDelayVector = InitialDelay*ones(1,num_grid_pts);
% %
% power = 30;
% power_converted = power*power_conversion;
% seq_stim_params.PowerVector = power_converted*ones(1,num_grid_pts);
% %2
% Duration = stim_duration;
% seq_stim_params.DurationVector = Duration*ones(1,num_grid_pts);
% % 
% numSpirals = 10;
% seq_stim_params.SpiralVector = numSpirals*ones(1,num_grid_pts);
% %
% Repetitions = 1;
% seq_stim_params.RepetitionsVector = Repetitions*ones(1,num_grid_pts);
% %Darcy sometimes recommends increasing 'Repetitions' and decreasing
% %Duration.  This changes the distribution of the spirals in time over the
% %cell.
% %
% InterPointDelay =  0.12;
% seq_stim_params.InterPointDelayVector = InterPointDelay*ones(1,num_grid_pts); 
% 
% % seq_stim_params
% %--------------------------------------------------------------------------
% createXmlFile_sequential_single_cell(xml_seq_path, seq_stim_params, stim_sequence);
% %--------------------------------------------------------------------------
% % Update prairie view repetitions based on num neurons to stim
% stim_time_per_neuron = InitialDelay/1000+InterPointDelay;
% len_seq_stim = numberNeurons*stim_time_per_neuron/60;
% 
% num_reps_seq_stim = ceil(expectedDuration*60*frameRate);
% 
% disp(['Number of Repetitions in PrairieView: ' num2str(num_reps_seq_stim)])
% disp('Update PrairieView to have at least this many reps'); 
% disp(['Length (min): ' num2str(expectedDuration)])
% 
% % if position of stim cells looks different, "smaller/bigger" check the
% % pixel size
% %--------------------------------------------------------------------------
% %DO: 
% % upload .gpl in MarkPoints (Top half)
% % upload .xml in MarkPoints (Bot half)
% % update T-series repetitions in Prairie View with above number
% % Make sure Voltage Recording has all channels enabled
% % Make sure you turn on the laser power and pmt's
% % while running paint the neurons
% % TODO: automate uploading
% %--------------------------------------------------------------------------
% 
% %%
% clear s
% % expt_str = 'gridstim'; %previously 'holostim' 
% expt_str = 'test_gridstim'; %previously 'holostim' 
% mask = roi_data.roi_mask;
% expectedLengthExperiment = ceil(num_reps_seq_stim*1.5); 
% HoloAcqnvsPrairie_v2(path_data, expt_str, mask, expectedLengthExperiment)


