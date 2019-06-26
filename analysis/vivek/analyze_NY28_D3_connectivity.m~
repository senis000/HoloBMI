
%%
%Code Paths
%--------------------------------------------------------------------------
addpath(genpath('/Users/vivekathalye/Dropbox/Code/holo_bmi')); 
addpath(genpath('/Users/vivekathalye/Dropbox/Code/analysis_util'))
addpath('/Users/vivekathalye/Dropbox/Code/export_fig/export_fig_11.19.15'); 
addpath('/Users/vivekathalye/Dropbox/Code/holobmi_git/HoloBMI/baseline_target_calibration')

%%
%Data Paths
%--------------------------------------------------------------------------
data_dir = '/Users/vivekathalye/Dropbox/Data/holo_may2019/NY28_190521'
exist(data_dir)

im_dir = '/Users/vivekathalye/Dropbox/Data/holo_may2019/NY28_190521/im/connectivity_pre/connectivity_pre_190521T223307-020'
exist(im_dir)

%%
%Tifs
tifList = dir(fullfile(im_dir, '*.tif'))
num_frames_saved = length(tifList); 

%%
%voltageRec
voltageRec_file = fullfile(im_dir, 'connectivity_pre_190521T223307-020_Cycle00001_VoltageRecording_001.csv'); 
exist(voltageRec_file)
row_offset = 2; 
voltageRec = csvread(voltageRec_file, 2, 0);
v_sampling_freq = 1000; %Hz, samples/sec
%%
%indexing voltage: stim time, frame time
%(2nd is power)
v_t_idx         = 1; 
v_laser_idx     = 8; 
v_frame_idx     = 3;

%Laser Idx: 8
%time stamp (ms): 1

%%
%Frame idxs, frame times
frame_detect_thresh = 2; %(max(v_frame)-min(v_frame))/2;

t = voltageRec(:,v_t_idx)/(v_sampling_freq*60); %minutes
v_frame = voltageRec(:,v_frame_idx);

frame_diff = conv(v_frame, [1 -1], 'valid');
% h = figure;
% hold on;
% plot(t, v_frame); 
% plot(t(2:end), frame_diff);
% title('Frame vs Frame Diff'); 

frame_idxs = find(frame_diff > frame_detect_thresh)+1; 
frame_times = t(frame_idxs);
frame_times_saved = frame_times(1:num_frames_saved); 
% vline(frame_times); 

%%
% Stim idxs, stim times
t = voltageRec(:,v_t_idx)/(v_sampling_freq*60); %minutes
stim_voltage = voltageRec(:,v_laser_idx);
% h = figure;
% plot(t, stim_voltage); 
% xlabel('time (min)')
% ylabel('stim'); 

%
%detect times of voltage: 
stim_detect_thresh = 0.03; %hand coded
stim_diff = conv(stim_voltage, [1 -1], 'valid');
h = figure;
hold on;
plot(t, stim_voltage);
% plot(t(2:end), stim_diff);
stim_idxs = find(stim_diff > stim_detect_thresh)+1;
stim_times  = t(stim_idxs);
vline(stim_times); 

%%
%Assign stim to the frame just after it occurs
num_stim_delivered =length(stim_times); 
stim_frame = ones(num_stim_delivered,1)*-1; 
for i = 1:num_stim_delivered
    tdiff_i = stim_times(i)-frame_times_saved; 
    %Want the frame just before stim time, so tdff_i > 0
    [Y,I] = min(abs(tdiff_i(tdiff_i > 0))); 
    stim_frame(i) = I(1);    
end

%%
%Which neurons were stimmed in Experiment:
target_info_cell = {fullfile(data_dir, 'BMI_target_info_20190521T224235.mat')}; 
tinfo = load(target_info_cell{1}); 

%Masks used to for the above cells: 
im_file_cell = {...
    fullfile(data_dir, 'red.mat')};
%Use: holoMaskRedGreen
%Load the red.mat, get the masks to extract the activity from the videos: 
load(im_file_cell{1}); 

%Plot colors 
plot_b = [0    0.4470    0.7410];
plot_o = [0.8500    0.3250    0.0980]; 
E_color = {plot_b, plot_o}; 
num_BMI = length(tinfo.E1_base) + length(tinfo.E2_base);

%%
%Find Order In Which Neurons were Stimmed
%xml used for MarkPoints
xml_file = fullfile(data_dir, 'connectivity_pre_holostim.xml'); 
% xml_file = fullfile(im_dir, 'connectivity_pre_190521T223307-020.xml'); 
exist(xml_file)

%Read the xml file as a string, look for point:
xmltext = fileread(xml_file);
stim_text_idxs = strfind(xmltext, 'Point ')

num_stim = length(stim_text_idxs);
stim_order = -1*ones(num_stim, 1); 
len_to_num = 6; %6 char after finding 'Point ' is the index of the neuron stimmed
for i = 1:num_stim
    stim_order(i) = str2num(xmltext(stim_text_idxs(i)+len_to_num));
%     xmltext(stim_text_idxs(i)+len_to_num)
end
%Only take the stims that were actually delivered
stim_order_delivered = stim_order(1:num_stim_delivered); 

stim_label = cell(num_stim_delivered,1); 
for i=1:num_stim_delivered
    stim_label{i} = num2str(stim_order_delivered(i)); 
end
%For Labeling which cell was stimmed

%To Do: 
%In future, save the order of stim in the BMI .mat for easier analysis.

%%
%Get mask for E2 stimmed neurons:
%holoMaskRedGreen from 'red.mat'
E2_mask     = zeros(size(holoMaskRedGreen)); 
for indn=1:length(tinfo.E2_base)
    auxmask = holoMaskRedGreen;
    auxmask(auxmask~=tinfo.E2_base(indn)) = 0;
    auxmask(auxmask~=0) = indn;
    E2_mask = auxmask + E2_mask;    
end
E2_strcmask = obtainStrcMaskfromMask(E2_mask);

%%
%Extract activity of stimmed neurons:
num_frames  = length(tifList); 
num_E2      = length(tinfo.E2_base); 
unitVals    = zeros(num_E2,num_frames);
tic
for frame_i=1:num_frames
    data_frame      = imread(fullfile(im_dir, tifList(frame_i).name));
    unitVals(:, frame_i)  = obtainRoi(data_frame, E2_strcmask, 1:num_E2);
end
toc

%%
%Plot activity and stim
t_as_frame_num = 1; 

n = unitVals.';
if(t_as_frame_num)
    t = 1:length(n); %frame_times(1:length(n)); 
else
    t = frame_times_saved; 
end
%(1:length(n))*0.03/60;
E_id = 2*ones(num_E2, 1);
[h, offset_vec] = plot_E_activity(t, n, E_id, E_color, 0);
hold on; 
if(t_as_frame_num)
    vline(stim_frame,  'r:', stim_label); 
else
    vline(stim_times,  'r:', stim_label); 
%     vline(stim_times, repmat({'r:'}, num_stim_delivered, 1), stim_label); 
end

%%
%Annotate possibilities: 
%fade non-stim neurons
%draw a faded contour
%light up the contour at the right moments
% stim box in a corner

fade_coeff = 0.1; 



% E2_mask_bin = E2_mask > 0; 
mask_bin = (E2_mask==2); 
mask_contrast = mask_bin; 
mask_contrast(mask_contrast) = 0.1; 

%%
% h = figure;
% imagesc(mask_bin); 
% colormap('gray'); 
% axis square

%%
%Load the entire video, and multiply with mask_contrast



%%
vid_sel = 1150:1250; 
num_sel = length(vid_sel); 
num_row = 512;
num_col = 512; 
A = zeros(num_row, num_col, num_sel);  
A_annotate = A; 
for i=1:num_sel
    frame_i                 = vid_sel(i); 
    data_frame              = imread(fullfile(im_dir, tifList(frame_i).name));
	A(:,:,i)          = data_frame; 
    A_annotate(:,:,i) = double(data_frame).*(mask_bin); 
%     A_annotate(:,:,frame_i) = data_frame.*mask_contrast; 
end

%%
%

%%
%frames to look at:
%1200, 1208

frame_i = 1200; 
%
clim = [0 6400]; 
%
i = find(vid_sel==frame_i); 
data_frame      = A(:,:,i); 
h = figure;
imagesc(data_frame, clim); 
title(['full ' num2str(frame_i)]); 
colormap('gray')
axis square

%
i = find(vid_sel==frame_i); 
data_frame      = A_annotate(:,:,i); 
h = figure;
imagesc(data_frame, clim);
title(['roi ' num2str(frame_i)]); 
colormap('gray')
axis square

%%
h = figure;
imagesc(data_frame)

%--------------------------------------------------------------------------
%%
h = figure; 
imagesc(E2_mask);

%%
h = figure;
hold on;
plot(unitVals(1,:)); 

%%
%unitVals: 4 x 1500
% [h, offset_vec] = plot_E_activity(t,n, E_id, E_color, offset)

%%


%%


%%
% h = figure;
% imagesc(E2_mask==3)
% axis square

%%
h = figure;
%Plot 



%%
num_stim_exec           = length(stim_times); 
stim_order_exec         = stim_order(1:num_stim_exec);  

% t = frame_times(1:length(n)); 
t = 1:length(n); 
for i = 1%:num_E2
    i_sel = (stim_order_exec == i);
%     stim_i = stim_times(i_sel); 
    stim_i = stim_frame(i_sel); 
    
    h = figure;
    plot(t, n(:,i)); 
    hold on;
    vline(stim_i); 
end


%%
num_stim_exec           = length(stim_times); 
stim_order_exec         = stim_order(1:num_stim_exec);  
for i = 4%:num_E2
    i_sel = (stim_order_exec == i);
    stim_i = stim_times(i_sel); 
    
    h = figure;
    plot(t, n(:,i)); 
    hold on;
    vline(stim_i); 
end

%%



%%
[h, offset_vec] = plot_E_activity(t, n, E_id, E_color)

%%
num_rows = 512;
num_cols = 512; 
num_frames = length(tifList); 
A = zeros(num_rows, num_cols, num_frames); 



%%
h = figure;
imagesc(A(:,:,288))
axis square