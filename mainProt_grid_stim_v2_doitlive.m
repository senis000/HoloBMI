
%%
%--------------------------------------------------------------------------
%BEFORE ANIMAL IN BOX:
%DO:
% Hook up BNCs: 
% 1) BMI solenoid, AI5
% 2) Monaco Trig, AI6
% 3) Frame Trig, AI7
% 4) Holo Trig PFI1

%%
% %Load example roi_data information for reference: 
% roi_data_path = fullfile('Z:\Vivek\holoBMI_2round\191119\NVI20\D15', 'roi_data.mat'); 
% d = load(roi_data_path)
% workspace_path = fullfile('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-13', 'workspace.mat'); 
% w = load(workspace_path)

%% DEFINE PATHS

cd G:\VivekNuria\Code\HoloBMI
%DEFINE PATH_DATA: 

% define Animal, day and folder where to save
%Animals: NVI20, NVI21
animal = 'NVI20'; day = 'D0';
folder = 'E:\holobmi_E\200214';
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
% REMEMBER TO TURN OFF PHASE OFFSET
% TURN OFF THE MANIPULATOR 
% TURN OFF AUTOSCALE
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
%--------------------------------------------------------------------------
%STIM RELATED CODE HERE
%--------------------------------------------------------------------------

%Define target location of stim
%Default to center of image: 
%Can manually select x, y based on roi_data
target.x = 258 %dim/2
target.y = 266 %dim/2
target.r = 10  %not super necessary for stim grid points.  needed for mark points generation

%%
%1) define grid in pixel space
% micronsPerPixel.x = .7266
% micronsPerPixel.y = .7266
%USER INPUT: 
dim = 512
x_extent_pixels = [-dim dim]/4
y_extent_pixels = [-dim dim]/4


x_extent_centered = target.x + x_extent_pixels
x_extent_centered(x_extent_centered < 1) = 1;
x_extent_centered(x_extent_centered > dim) = dim;

y_extent_centered = target.y + y_extent_pixels;
y_extent_centered(y_extent_centered < 1) = 1;
y_extent_centered(y_extent_centered > dim) = dim;

step_microns = 10
x_step = round(step_microns/micronsPerPixel.x);
y_step = round(step_microns/micronsPerPixel.y);
%
x_step_vec = x_step:x_step:x_extent_pixels(2);
y_step_vec = y_step:y_step:y_extent_pixels(2);
%

remove_x = find(...
    (target.x+x_step_vec) < x_extent_centered(1) | ...
    (target.x+x_step_vec) > x_extent_centered(2))
x_step_vec(remove_x) = [];

remove_y = find(...
    (target.y+y_step_vec) < y_extent_centered(1) | ...
    (target.y+y_step_vec) > y_extent_centered(2))
y_step_vec(remove_y) = [];

x_mesh_flat = target.x + [0 -x_step_vec zeros(size(y_step_vec))]; 
y_mesh_flat = target.y + [0 zeros(size(x_step_vec)) -y_step_vec];
 
ctr_idx = find(x_mesh_flat == target.x & y_mesh_flat == target.y); 
num_grid_pts = length(x_mesh_flat) 


%PLOTTING:
% h = figure;
% plot(x_grid, '.-')

% h = figure; 
% imagesc(X_mesh)
% colorbar
% axis square
% h = figure; 
% imagesc(Y_mesh)
% colorbar
% axis square

h = figure;
hold on
scatter(x_mesh_flat, y_mesh_flat, 'x'); 
scatter(x_mesh_flat(ctr_idx), y_mesh_flat(ctr_idx), 'r', 'x'); 
xlim([1 dim])
ylim([1 dim])
axis square
%Reminder, mapping from x,y to image space: 
%x = column, y = row

%%
%--------------------------------------------------------------------------
%2) map grid to stim points (GPL file): 
%DEFINE GPL FILE
spiral_size_conversion = 1/49; 
%a coefficient needed to accurately load desired ROI size
%Empirically measured
initSpiralSize_um = 14; 
initSpiralSize = spiral_size_conversion*initSpiralSize_um;

stim_duration       = 5;  %(ms)

init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', stim_duration, ...
    'SpiralSize', initSpiralSize, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals

markpoints_data = repmat(init_markpoints, [num_grid_pts 1]); 
createGplFile_v2(savePath, markpoints_data, x_mesh_flat, y_mesh_flat, posz, NaN, px, zoom);

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
num_grid_stim_per_target_stim = 10;

num_stims_expected = num_reps*(num_grid_pts + ceil(num_grid_pts/num_grid_stim_per_target_stim))

rand_order = randperm(num_grid_pts); 
stim_sequence = [ctr_idx ctr_idx]

for i = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
    sel_from_rand_order = ...
        (1:num_grid_stim_per_target_stim)+(i-1)*num_grid_stim_per_target_stim;   
    sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
    grid_idxs = rand_order(sel_from_rand_order); 
    for j = 1:length(grid_idxs)
        stim_sequence = [stim_sequence grid_idxs(j) grid_idxs(j)];
    end
    stim_sequence = [stim_sequence ctr_idx ctr_idx];
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
%Confirmation:
h= figure;
hist(stim_sequence, 1:num_grid_pts)
title(['num stims per grid point, ctr idx: ' num2str(ctr_idx)])
disp('center point:')


%% XML: Sequential Single Cell Stim

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
time_between_stims  = 10000; %(ms)
expectedDuration = (stim_duration + time_between_stims)*num_stims/(60*1000) %min 

%%
xml_seq_path = fullfile(savePath, 'grid_stim.xml'); 
power_conversion = 0.004; %0.2 -> 50, 0.4->100

seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
%
InitialDelay = time_between_stims; %(ms) time bw stim delivery
seq_stim_params.InitialDelayVector = InitialDelay*ones(1,num_stims);
%
power = 20;
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
%ToDO: stim analysis of the results

%%
%For fast check
%Just up and right, 
%
