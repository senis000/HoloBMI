function count_intra_block_event(block_data, block_event, param)
% -first half, last half # events
% -cumulative # of hits over time
% -hit rate with sliding window
% input:
% -param:
% -- frame_rate (frames/s)
% -- save_bool
% -- save_dir
% -- sliding_win (in minutes)
% -- cumsum_times
% -- win_times : evaluate sliding window at these times
% -- num_equal_paritions

r = struct(...
    'num_events', [], ...
    'event_rate', [], ...
    'num_events_first_last', [], ...
    'cum_cnt', [], ...
    't_cum_cnt', [], ...
    'event_rate_sliding', [], ...
    't_event_rate_sliding', []); 

%--
time_scale = 1/(frame_rate*60); %minutes per frame
block_len = size(block_data.dff, 2)*time_scale; 

%--
%num_events, event_rate
r.num_events = length(block_event.data); 
r.event_rate = r.num_events/block_len;

%--
%num_events_first_last
block_half_time = block_len/2;
r.num_events_first_last = ...
    [sum(block_event.data <= block_half_time) ...
    sum(block_event.data >= block_half_time)];

%--
%cum_cnt, t_cum_cnt
%To Do: function to interpolate timeseries t1,y1 to t2,y2 using zero order hold
r.cum_cnt = 1:r.num_events;
r.t_cum_cnt = block_event.data;

%--
%hit_rate_sliding
%For loop over win_partitions:
r_win = zeros(param.win_times,1); 
for i = 1:length(param.win_times)
    stop_i = param.win_times(i); 
    start_i = stop_i - param.sliding_iwn; 
    r_win(i) = sum(block_event.data>=start_i & block_event.data<=stop_i);
    r_win(i) = r_win(i)/(stop_i-start_i); 
end
r.event_rate_sliding    = r_win;
r.t_event_rate_sliding  = ;


% num_files   = size(block_event,1);
% num_events  = size(block_event,2);
% 
% % result = repmat(struct(...
% %     'num_events', [], ...
% %     'hit_rate', [], ...
% %     'num_events_first_last', [], ...
% %     'cum_perc', [], ...
% %     'hit_rate_sliding', [] ...
% %     ), [num_files num_events]); 
% 
% 
% for i_file = 1:num_files
%     for i_event = 1:num_events
%         block_i = block_data(i_file); 
%         event_i = block_event(i_event); 
%         result(file_i, event_i) = count_event(block_i, event_i); 
%     end
% end

end

% function count_event(block_data, event_data)
% 
% 
% end

