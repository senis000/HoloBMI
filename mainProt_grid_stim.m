
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
%1) define grid in pixel space
% micronsPerPixel.x = .7266
% micronsPerPixel.y = .7266

dim = 512

target.x = dim/2
target.y = dim/2
target.r = 10  %not super necessary for stim grid points.  needed for mark points generation

x_extent_pixels = [-dim/8 dim/8]
y_extent_pixels = [-dim/8 dim/8]
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
x_grid = target.x + [sort(-x_step_vec) 0 x_step_vec]; 
remove_x = find(...
    x_grid < x_extent_centered(1) | ...
    x_grid > x_extent_centered(2))
x_grid(remove_x) = []

y_grid = target.y + [sort(-y_step_vec) 0 y_step_vec]; 
remove_y = find(...
    y_grid < y_extent_centered(1) | ...
    y_grid > x_extent_centered(2))
y_grid(remove_y) = []
%Sort is just used for ease of plotting the results
%
[x_mesh, y_mesh] = meshgrid(x_grid, y_grid);
%
x_mesh_flat = x_mesh(:); 
y_mesh_flat = y_mesh(:); 
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
initSpiralSize_um = 16; 
initSpiralSize = spiral_size_conversion*initSpiralSize_um/micronsPerPixel.x;

init_markpoints = struct(...
    'UncagingLaserPower', 0.4, ...
    'Duration', 100, ...
    'SpiralSize', initSpiralSize, ...
    'SpiralRevolutions', 10); 
%Darcy recommends 5-10 spirals

markpoints_data = repmat(init_markpoints, [num_grid_pts 1]); 
createGplFile_v2(savePath, markpoints_data, x_mesh_flat, y_mesh_flat, posz, roi_data.r, px, zoom);

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
stim_sequence = []
for i = 1:num_reps
    rand_order = randperm(num_grid_pts);
    for j = 1:ceil(num_grid_pts/num_grid_stim_per_target_stim)
        sel_from_rand_order = ...
            (1:num_grid_stim_per_target_stim)+(j-1)*num_grid_stim_per_target_stim;
        %Remove invalid idxs
        sel_from_rand_order(sel_from_rand_order > num_grid_pts) = []; 
        grid_idxs = rand_order(sel_from_rand_order); 
        stim_sequence = [stim_sequence ctr_idx grid_idxs]
    end
end

%Confirmation:
h= figure;
hist(stim_sequence, 1:num_grid_pts)
title(['num stims per grid point, ctr: ' num2str(ctr_idx)])
disp('center point:')

%% XML: Sequential Single Cell Stim
xml_seq_path = fullfile(savePath, 'grid_stim.xml'); 
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
time_between_stims = 10000; %(ms)



power_conversion = 0.004; %0.2 -> 50, 0.4->100

seq_stim_params.UncagingLaser = "Monaco"; 
seq_stim_params.AllPointsAtOnce = "False"
seq_stim_params.Iter = 1; %how many times to go through and stim each cell.
seq_stim_params.IterDelay = 1000; %Time (ms) between iterations
%
InitialDelay = time_between_stims; %(ms) time bw stim delivery
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

