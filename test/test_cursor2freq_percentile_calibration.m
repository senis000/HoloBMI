
%Test mapping cursor values to frequencies: 

%%
data = load(cal.paths.cal_all)
%%
h = figure;
hist(data.cursor_obs, 50)

%%

cal.fb.cursor_min
cal.fb.cursor_max

%%
num_intervals = 7
cursor_equal_interval_boundaries = ...
    linspace(cal.fb.cursor_min, cal.fb.cursor_max, num_intervals+1); 
freq_equal_interval_boundaries = ...
    cursor2audio_freq(cursor_equal_interval_boundaries, cal)

%%
h = figure;
plot(freq_equal_interval_boundaries, '.-', 'MarkerSize', 15)


%%
cursor_interval_boundaries 
cursor_interval_boundaries = [cal.fb.cursor_min; cal.fb.cursor_max]

%%
%1) calculate percentile vector: 
num_intervals = 2; 
num_percentiles = num_intervals - 1; 
percentiles_with_ends = linspace(0,100, num_percentiles+2)

% h = figure;
% vline(percentiles_with_ends); 

percentile_without_ends = percentiles_with_ends(2:end-1); 

% h = figure;
% vline(percentile_without_ends); 

% cursor_boundaries = 

%%
test = prctile(data.cursor_obs, percentiles_with_ends)
h = figure;
hold on;
hist(data.cursor_obs, 50); 
vline(test); 

%%
prc_val_without_ends    = prctile(data.cursor_obs, percentile_without_ends)
prc_val_with_ends       = [cal.fb.cursor_min prc_val_without_ends cal.fb.cursor_max]
h = figure;
hold on;
hist(data.cursor_obs, 50); 
vline(prc_val_with_ends); 

%%
%option 1: make it general for a desired tone distribution
%option 2: just split it in 2
%for this, just calculate median, calculate two exponentials. 

%%
fake_cursor = linspace(cal.fb.cursor_min, cal.fb.cursor_max, 10000);

%%
fake_fb = cursor2audio_freq_middle_match(fake_cursor, cal);

h = figure;
plot(fake_cursor, fake_fb); 

%%
cal_test = cal; 
cal_test.fb.target_low_freq = 0; 
%%
[fake_fb, fb_combine, bin_sel]  = cursor2audio_freq_middle_match(fake_cursor, cal);

%%
ctrl_fb = cursor2audio_freq(fake_cursor, cal);

h = figure;
% plot(fake_cursor, ctrl_fb); 
% hold on; 
plot(fake_cursor, fake_fb); 

%%
h = figure;
hold on;
% plot(fake_cursor, fb_combine(:,1)); 
plot(fake_cursor, fb_combine(:,2)); 

%%
test = zeros(size(fb_combine,1), 1); 
for i =1:length(test)
    if fake_cursor(i) <= cal.fb.cursor_middle
        test(i) = fb_combine(i,2); 
    else
        test(i) = fb_combine(i,1); 
        
    end
end

%%
bin_sel_test = (-fake_cursor >= cal.fb.cursor_middle) + 1

%%
fake_fb_test = zeros(size(fb_combine, 1), 1); 
for i =1:length(fb_combine)
    fake_fb_test(i) = fb_combine(i,bin_sel_test(i)); 
end
h = figure; 
plot(fake_fb_test); 

%%
h = figure;
hold on;
plot(fake_cursor, test); 
plot(fake_cursor, fake_fb); 

%%
% freq_middle: 8.2608e+03
% freq_middle: 1.0677e+04
%%
h = figure;
plot(fake_cursor, ctrl_fb); 
