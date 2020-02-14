%test_load_roi_data

%%
roi_data_path = fullfile('Z:\Vivek\holoBMI_2round\191119\NVI20\D15', 'roi_data.mat'); 
d = load(roi_data_path)

%%
workspace_path = fullfile('Z:\Ines\2photon\Training\ISO+BMI 1\NY127\2019-12-13', 'workspace.mat'); 
w = load(workspace_path)

%%
h = figure;
imshow(d.roi_mask)

%%
sel_idx = 1
roi_im = d.roi_data.roi_bin_cell{sel_idx};
h = figure;
imshow(roi_im); 

x = d.roi_data.x(sel_idx)
y = d.roi_data.y(sel_idx)
%%