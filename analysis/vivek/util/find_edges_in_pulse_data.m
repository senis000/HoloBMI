function [edges] = find_edges_in_pulse_data(data, edge_thresh_coeff)
%rising idx: the idx before the rise happens
%falling idx: the idx before the fall happens
[data_norm, data_min, data_max, data_range] = range_norm(data, 0, 100);
edge_thresh = edge_thresh_coeff; 
data_peaky = data; 
data_peaky(data_peaky<edge_thresh) = 0; 

%Finite difference
diff_filt = [1 -1];
data_diff = conv(data_peaky, diff_filt, 'same'); 

rise_bool = data_diff;
rise_bool(data_diff <= edge_thresh)     = 0;
rise_bool(data_diff > edge_thresh)      = 1;

fall_bool = data_diff;
fall_bool(data_diff >= -edge_thresh)    = 0;
fall_bool(data_diff < -edge_thresh)     = 1;

edges.rise = find(rise_bool); 
edges.fall = find(fall_bool);


% def find_edges_in_pulse_data(data, edge_thresh_coeff):
% 	#data: 
% 		#assumes pulse data ranges from 0 to 'amp'
% 	#an edge is detected if data changes by >= amp*edge_thresh_coeff
% 	#edge_thresh_coeff:
% 		# value from 0 to 1
% 
% 	edges = {}
% 
% 	amp = np.amax(data)
% 	edge_thresh = amp*edge_thresh_coeff
% 
% 	data_peaky = np.copy(data)
% 	data_peaky[data<edge_thresh] = 0
% 
% 	#Finite difference 
% 	diff_filt = np.array([1, -1])
% 	data_diff = np.convolve(data_peaky, diff_filt)
% 
% 	#Rise:
% 	rise_bool = np.copy(data_diff)
% 	rise_bool[data_diff <= edge_thresh] = 0
% 	rise_bool[data_diff > edge_thresh] = 1
% 	edges['rise'] = np.where(rise_bool)[0]
% 
% 	#Fall: 
% 	fall_bool = np.copy(data_diff)
% 	fall_bool[data_diff >= -edge_thresh] = 0
% 	fall_bool[data_diff < -edge_thresh] = 1
% 	edges['fall'] = np.where(fall_bool)[0]
% 
% 	return edges