
%%
im_home = 'D:\DATA\grid_stim\200311\NVI20\D0\im'

%%
file = 'fine_gridstim-004'

data_dir = fullfile(im_home, file); 
v_file = [file '_Cycle00001_VoltageRecording_001.csv']
v_path = fullfile(data_dir, v_file)
exist(v_path)
tic
voltageRec = readmatrix(v_path);
toc
disp('loaded voltageRec!'); 

frame_chan = 1+1+1;
stim_chan = 1+6+1; 

plot_bool = 1
if plot_bool
    h = figure;
    plot(voltageRec(:, frame_chan)); 
    title('frame pulse'); 

    h = figure;
    plot(voltageRec(:,stim_chan))
    title('stim pulse'); 
end

min_percentile = 0; 
max_percentile = 100; 
[frame_pulse, frame_min, frame_max, frame_range] = range_norm(voltageRec(:, frame_chan), max_percentile, min_percentile);

min_percentile = 0; 
max_percentile = 100; 
[stim_pulse, stim_min, stim_max, stim_range] = range_norm(voltageRec(:, stim_chan), max_percentile, min_percentile);

edge_thresh_coeff = 0.5; 
stim_edge   = find_edges_in_pulse_data(stim_pulse, edge_thresh_coeff);
frame_edge  = find_edges_in_pulse_data(frame_pulse, edge_thresh_coeff);
% [edges] = find_edges_in_pulse_data(data, edge_thresh_coeff)

%
%Map the stim_idx to the frame_idx that precedes 
stim_frame = []; %zeros(size(stim_edge.rise));
for i=stim_edge.rise(:)'
    i_diff = i-frame_edge.rise; 
    i_diff(i_diff > 0) = -1e7; 
    [M,I] = max(i_diff);
    stim_frame = [stim_frame I]; 
end
%I visually verified that the target cell responds to the stimulation in
%imagej

disp('RESULT:')
file
size(frame_edge.rise)

%%
%fine_gridstim-000: 
%22013.  
%lastGoodFrame=10
%
%fine_gridstim-001: 
%22012
%lastGoodFrame=11
%
%fine_gridstim-002:
%22004
%lastGoodFrame=0
%
%fine_gridstim-003:
%22013
%lastGoodFrame=10
%
%fine_gridstim-004:
%22012
%lastGoodFrame=12

