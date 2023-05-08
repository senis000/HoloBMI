function result = grid_stim_files2psth_data(psth_win, tif_dir, grid_data_path, bg_im_path, v_path, save_dir)

load(grid_data_path) %contains data about the grid for stimulation
%x_mesh_flat, y_mesh_flat have the idxs for grid points
%ctr_idx: idx for the ctr point

im_bg = imread(im_bg_path); 
