data_dir = 'G:\vivek\DATA'
data_path = fullfile(data_dir, 'playback_base111619.mat')

%%
test = load(data_path)

%%
h = figure;
hold on;
hist(test.cursor_obs, 50)
vline(test.cal.target.E2_hit_cal.T); 
vline(test.cal.target.E1_hit_cal.T); 

%%
h = figure;
hist(log(test.fb_obs), 50)