
%%
%Code paths to add: 
addpath('D:\Dropbox\Code\export_fig\export_fig_3.3.20'); 
%%
%file to load: 
clear all
date_str    = '200311'
file_analyze = 'fine_gridstim-004' %'stim-005'
%TODO: select the files to load
%fine_gridstim-000
%fine_gridstim-001
%fine_gridstim-002
%fine_gridstim-003
%fine_gridstim-004
%
%coarse_gridstim-000
%coarse_gridstim-001
%coarse_gridstim-002
%coarse_gridstim-003
%coarse_gridstim-004

%tif_dir     = ['D:\DATA\grid_stim\' date_str '\NVI20\D0\im\' file_analyze] %coarse_grid-000'
tif_dir     = 'D:\DATA\grid_stim\200311\NVI20\D0\im\fine_gridstim-004'
data_dir    = 'D:\DATA\grid_stim\200311\NVI20\D0'
load_path   = fullfile(data_dir, 'fine_grid.mat')
im_bg_path          = fullfile('D:\DATA\grid_stim\200311\NVI20\D0\redgreen', 'green.tif'); 
voltageRec_path     = ...
    fullfile('D:\DATA\grid_stim\200311\NVI20\D0\im\fine_gridstim-004', 'fine_gridstim-004_Cycle00001_VoltageRecording_001.csv');
%fullfile('D:\DATA\grid_stim\200311\NVI20\D0\im\fine_gridstim-000', 'fine_gridstim-000_Cycle00001_VoltageRecording_001.csv');  
%fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']); 

%['D:\DATA\grid_stim\' date_str '\NVI20\D0\coarse_grid_data' %['D:\DATA\grid_stim\200214\NVI20\D0\data_' file_analyze];

save_dir    = fullfile('D:\Dropbox\Data\grid_stim', date_str, file_analyze); 
mkdir(save_dir); 




%
screensize = get(0, 'ScreenSize')

%clear
load(load_path)

%%
% animal              = 'NVI20'
% date_str            = '200312_H'
% file_analyze        = 'fine_gridstim-001' %'coarse_gridstim-004' %'fine_gridstim-004'
% psth_win_samples    = [-8 8]*30;
% stim_lag_frames     = 0;
% data_dir            = fullfile('D:\DATA\grid_stim', date_str, animal, 'D0');
% exist(data_dir)
% grid_data_path      = fullfile(data_dir, 'fine_grid.mat'); 
% % grid_data_path      = fullfile(data_dir, 'coarse_grid.mat'); 
% exist(grid_data_path)
% im_bg_path          = fullfile(data_dir, 'redgreen', 'green.tif');
% exist(im_bg_path)
% tif_dir             = fullfile(data_dir, 'im', file_analyze);
% exist(tif_dir)
% v_path              = fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']);
% exist(v_path)
% frame_chan          = 3; 
% stim_chan           = 8; 
% save_dir            = fullfile('D:\Dropbox\Data\grid_stim', date_str, file_analyze); 
% mkdir(save_dir); 
% save_path           = fullfile(save_dir, 'psth_data.mat');

%%
data_dir            = fullfile('D:\DATA\grid_stim', date_str, animal, 'D0');
exist(data_dir)
%%
grid_data_path      = fullfile(data_dir, 'coarse_grid.mat'); 
coarse_test           = load(grid_data_path)

%%
grid_data_path      = fullfile(data_dir, 'fine_grid.mat'); 
fine_test           = load(grid_data_path)

%%
%Let's do this: 
%200310
animal              = 'NVI20'
date_str            = '200311'
file_analyze        = 'coarse_gridstim-003' %'coarse_gridstim-004' %'fine_gridstim-004'
psth_win_samples    = [-8 8]*30;
stim_lag_frames     = 0;
data_dir            = fullfile('D:\DATA\grid_stim', date_str, animal, 'D0');
exist(data_dir)
% grid_data_path      = fullfile(data_dir, 'fine_grid.mat'); 
grid_data_path      = fullfile(data_dir, 'coarse_grid.mat'); 
exist(grid_data_path)
im_bg_path          = fullfile(data_dir, 'redgreen', 'green.tif');
exist(im_bg_path)
tif_dir             = fullfile(data_dir, 'im', file_analyze);
exist(tif_dir)
v_path              = fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']);
exist(v_path)
frame_chan          = 3; 
stim_chan           = 8; 
save_dir            = fullfile('D:\Dropbox\Data\grid_stim', date_str, file_analyze); 
mkdir(save_dir); 
save_path           = fullfile(save_dir, 'psth_data.mat');

%
grid_stim_files2psth_data(psth_win_samples, stim_lag_frames, tif_dir, grid_data_path, im_bg_path, v_path, frame_chan, stim_chan, save_path)
save_bool = 1; 
plot_grid_psth_data(save_path, save_dir, save_bool);


%%
%Let's do this: 
%200310

animal              = 'NVI20'
date_str            = '200310_finegrid'
file_analyze        = 'fine_gridstim-002' %'coarse_gridstim-004' %'fine_gridstim-004'
psth_win_samples    = [-8 8]*30;
stim_lag_frames     = 0;
data_dir            = fullfile('D:\DATA\grid_stim', date_str, animal, 'D0');
exist(data_dir)
grid_data_path      = fullfile(data_dir, 'fine_grid.mat'); 
% grid_data_path      = fullfile(data_dir, 'coarse_grid.mat'); 
exist(grid_data_path)
im_bg_path          = fullfile(data_dir, 'redgreen', 'green.tif');
exist(im_bg_path)
tif_dir             = fullfile(data_dir, 'im', file_analyze);
exist(tif_dir)
v_path              = fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']);
exist(v_path)
frame_chan          = 3; 
stim_chan           = 8; 
save_dir            = fullfile('D:\Dropbox\Data\grid_stim', date_str, file_analyze); 
mkdir(save_dir); 
save_path           = fullfile(save_dir, 'psth_data.mat');

%
grid_stim_files2psth_data(psth_win_samples, stim_lag_frames, tif_dir, grid_data_path, im_bg_path, v_path, frame_chan, stim_chan, save_path)
save_bool = 1; 
plot_grid_psth_data(save_path, save_dir, save_bool);


%%
%Let's do this: 
%200310

animal              = 'NVI20'
date_str            = '200310_coarsegrid'
file_analyze        = 'coarse_gridstim-001' %'coarse_gridstim-004' %'fine_gridstim-004'
psth_win_samples    = [-8 8]*30;
stim_lag_frames     = 0;
data_dir            = fullfile('D:\DATA\grid_stim', date_str, animal, 'D0');
exist(data_dir)
% grid_data_path      = fullfile(data_dir, 'fine_grid.mat'); 
grid_data_path      = fullfile(data_dir, 'coarse_grid.mat'); 
exist(grid_data_path)
im_bg_path          = fullfile(data_dir, 'redgreen', 'green.tif');
exist(im_bg_path)
tif_dir             = fullfile(data_dir, 'im', file_analyze);
exist(tif_dir)
v_path              = fullfile(tif_dir, [file_analyze '_Cycle00001_VoltageRecording_001.csv']);
exist(v_path)
frame_chan          = 3; 
stim_chan           = 8; 
save_dir            = fullfile('D:\Dropbox\Data\grid_stim', date_str, file_analyze); 
mkdir(save_dir); 
save_path           = fullfile(save_dir, 'psth_data.mat');

%
grid_stim_files2psth_data(psth_win_samples, stim_lag_frames, tif_dir, grid_data_path, im_bg_path, v_path, frame_chan, stim_chan, save_path)
save_bool = 1; 
plot_grid_psth_data(save_path, save_dir, save_bool);


%%
test_dir = fullfile(save_dir, 'test'); 
mkdir(test_dir); 
save_bool = 0; 
plot_grid_psth_data(save_path, test_dir, save_bool)










%%
test = load(save_path)

%%


test_dir
% plot_grid_psth_data(save_path, test_dir)
