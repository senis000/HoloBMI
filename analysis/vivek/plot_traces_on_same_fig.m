function [min_val, max_val] = ...
    plot_traces_on_same_fig(x_data, traces, err_bool, traces_err, trace_color)
%7.20.18
%traces: num_neurons x num_samples
%traces_err: num_neurons x num_samples: (for plotting sem)
%trace_color: cell of strings indicating plot color for the trace

vert_offset_over_neuron = -(cumsum(max(traces, [], 2) + 5));
max_val = 0; 
min_val = -max(vert_offset_over_neuron); 
num_var = size(traces,1); 

% h = figure; 
% hold on;
for i = 1:num_var
    plot_data = traces(i,:)+vert_offset_over_neuron(i);
    plot(x_data, plot_data, trace_color{i});
    min_val = min(min_val, min(plot_data)); 
    if(err_bool)
        errbar(x_data, plot_data, traces_err(i,:), trace_color{i});
        min_val = min(min_val, min(plot_data-traces_err(i,:)));
    end
     
end
