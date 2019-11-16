function [freq] = cursor2audio_freq_middle_match(cursor, cal)
%Cursor is the neural control signal.  From old work, this would be: 
%sumE2-sumE1
%
%Auditory feedback calculation
%freq = a*exp(b*(cursor_trunc-cursor_min))
% freq_min = a =  a*exp(0)
% freq_max = a*exp(b*(cursor_max-cursor_min))
% b = log(freq_max/a)/(cursor_max-cursor_min)
%param: cursor_min, cursor_max
%%
%If params.fb.target_low_freq == 1, then we negate cursor and cursor_min,
%cursor_max
%
%This is because cursor is E2-E1.  Our target is gonna be positive.  

%%
% %Debug:
% cal.fb.target_low_freq = 1; 
% cal.fb.freq_max = 19000; 
% cal.fb.freq_min = 6000; 
% 
% cal.fb.cursor_min    = -1;
% cal.fb.cursor_max    = 1;
% cal.fb.cursor_range  = 2; 
% 
% cal.fb.a = cal.freq_min;
% cal.fb.b = (log(cal.freq_max) - log(cal.a))/cal.cursor_range; 
% 
% cursor = -4:0.1:2; 

%%
%Handle target -> freq:

%%

%Approach: Calculate both cursors and select the right one.


if cal.fb.target_low_freq == 1
    %This means cursor up makes auditory freq go down:
    bin_sel         = (cursor <= cal.fb.cursor_middle)+1; 
    cursor          = -cursor; 
    cursor_min      = -cal.fb.cursor_max;
    cursor_middle   = -cal.fb.cursor_middle; 
    cursor_max      = -cal.fb.cursor_min;
    
else
    %This means cursor up makes auditory freq go up:
    bin_sel         = (cursor >= cal.fb.cursor_middle)+1; 
    cursor_min      = cal.fb.cursor_min;
    cursor_middle   = cal.fb.cursor_middle; 
    cursor_max      = cal.fb.cursor_max;
end



%%
%lower: cursor_min -> cursor_middle
%upper: cursor_middle -> cursor_max

bot_lower           = cursor_min; 
top_lower           = cursor_middle; 
cursor_trunc_lower  = max(cursor, bot_lower); 
cursor_trunc_lower  = min(cursor_trunc_lower, top_lower); 
cursor_trunc_lower  = cursor_trunc_lower - bot_lower; 

bot_upper           = cursor_middle; 
top_upper           = cursor_max; 
cursor_trunc_upper  = max(cursor, bot_upper); 
cursor_trunc_upper  = min(cursor_trunc_upper, top_upper); 
cursor_trunc_upper  = cursor_trunc_upper - bot_upper; 

% cursor_trunc    = max(cursor, cursor_min); 
% cursor_trunc    = min(cursor_trunc, cursor_max); 
% cursor_trunc    = cursor_trunc - cursor_min; 

%%
freq_lower      = ...
    cal.fb.a_bin(1)*exp(cal.fb.b_bin(1)*(cursor_trunc_lower));
freq_lower      = freq_lower(:); 

freq_upper      = ...
    cal.fb.a_bin(2)*exp(cal.fb.b_bin(2)*(cursor_trunc_upper));
freq_upper      = freq_upper(:); 

freq_combine    = [freq_lower freq_upper]; 

% freq = zeros(size(freq_combine, 1), 1); 
% for i =1:length(freq_combine)
%     freq(i) = freq_combine(i,bin_sel(i)); 
% end
I = (1:size(freq_combine,1)).'; 
J = reshape(bin_sel, [], 1); 
k = sub2ind(size(freq_combine), I,J); 
freq = freq_combine(k); 
freq = double(freq); 

