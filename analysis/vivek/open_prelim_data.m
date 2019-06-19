%%
%once we understand the data, we can analyze neural dynamics from baseline
%and post
%mean, covariance, dynamics (linear)

data_dir = '/Users/vivekathalye/Dropbox/Data/holo_bmi_prelim'
data_path = fullfile(data_dir, 'vivek_prelim.mat')
load(data_path)

%%
plot_dir = fullfile('/Users/vivekathalye/Dropbox/Data/holo_bmi_prelim', 'plots'); 
mkdir(plot_dir)

%%
% signals=neurons values
% Online frames = actual number of frames (during online analysis it may jump from frame 1 to frame 3 )
% datainfo: info about how we got the data  (also the index of ensemble neurons)
% neurons: neurons value
% online_data: value that was used for online-holoBMI
% reward: where vta stim and vta_real where it was "volitional" or vta_holo when it was due to holographic stim

%Vivek: why is the size of OnlineData less than in "neurons" ?  How do I
%plot the online data against the saved events?  
size(OnlineData)
size(neurons.F)
%%
% datainfo = 
% 
%                     FPS: 15.2000
%                      E1: [1 3 4 7]
%                      E2: [2 5 6 8]
%               BL_frames: [1x15567 double]
%             Post_frames: [1x15266 double]
%         HoloStim_frames: [1x154 double]
%             Holo_frames: [1x37832 double]
%     which_online_neuron: [8x1 double]

%%
% max(datainfo.BL_frames)
% 
% ans =
% 
%        15567
% 
% max(datainfo.HoloStim_frames)
% 
% ans =
% 
%        53146


%%
% indirect_idxs = ;

%%
disp('end of baseline:')
datainfo.BL_frames(end)

disp('end of holo:')
datainfo.Holo_frames(end)

disp('end of post:')
datainfo.Post_frames(end)

%%
disp('baseline length in frames: ')
datainfo.BL_frames(end)-datainfo.BL_frames(1)

disp('stim period frames: ')
datainfo.Holo_frames(end)-datainfo.Holo_frames(1)

disp('post stim period frames: ')
datainfo.Post_frames(end)-datainfo.Post_frames(1)

%%
disp('baseline length in time (min): ')
(datainfo.BL_frames(end)-datainfo.BL_frames(1))/(datainfo.FPS*60)

disp('stim period time: ')
(datainfo.Holo_frames(end)-datainfo.Holo_frames(1))/(datainfo.FPS*60)

disp('post stim period time: ')
(datainfo.Post_frames(end)-datainfo.Post_frames(1))/(datainfo.FPS*60)


% baseline length in time (min): 
% 
% ans =
% 
%    17.0680
% 
% stim period time: 
% 
% ans =
% 
%    41.4814
% 
% post stim period time: 
% 
% ans =
% 
%    16.7379

%%
num_holostim = length(datainfo.HoloStim_frames)
num_vta_holo = length(reward.VTAholo)
num_unrewarded_holostim = num_holostim - num_vta_holo


%%
num_frames1 = length(datainfo.BL_frames) + length(datainfo.Holo_frames) + length(datainfo.Post_frames)
num_frames2 = size(neurons.F, 2)

datainfo.Holo_frames(end)
reward.VTAreal

%%
sum((reward.VTAreal <= datainfo.BL_frames(end)))

sum((reward.VTAreal >= datainfo.Holo_frames(1)) & (reward.VTAreal <= datainfo.Holo_frames(end)))

sum(reward.VTAreal > datainfo.Holo_frames(end))

%%
% neurons = 
% 
%          F: [90x68665 single]
%     spikes: [90x68665 single]
%       Fraw: [90x68665 single]
%       Sraw: [90x68665 single]

%F = dff
%Fraw = raw fluorescence
%Sraw = from caiman
%spikes = estimate from caiman

%%
% reward = 
% 
%         VTA: [174x1 double]
%     VTAreal: [36x1 double]
%     VTAholo: [138x1 double]

%%
%1) Let us plot the neural data over time
%2) Mark Holo Events
%3) Mark VTA stim events

%%
%Step 1: 
num_chan = size(neurons.F,1)
E1_chan = datainfo.which_online_neuron(datainfo.E1)
E2_chan = datainfo.which_online_neuron(datainfo.E2)
E2_chan = E2_chan(2:end);

indirect_chan = 1:num_chan;
indirect_chan([E1_chan; E2_chan]) = [] 
%cut off the first E2 channel because there's a repeat

% E1_chan =
% 
%     74
%     61
%     26
%     23
% 
% 
% E2_chan =
% 
%     43
%     43
%     19
%     39
%
%From Nuria:
% Yes there was a mistake that day and we had twice the same neuron in E2

%%
% test = neurons.Sraw(E1_chan,:);
% size(test)

h = figure;
hold on
for i = 1:4
    plot(neurons.F(E1_chan(i),:))
end

%%
%Assign neural source: 
neural = neurons.Fraw;

num_frames = size(neural, 2); 
x_data = (1:num_frames);
frame2time = (1/(datainfo.FPS*60));
%Plot E1 and E2
num_E1 = length(E1_chan)
num_E2 = length(E2_chan)
num_indirect = length(indirect_chan); 

bmi_colors.E1 = 'r';
bmi_colors.E2 = 'k';
bmi_colors.indirect = 'm';

plot_x_in_time = 1; 
plot_holo_stim_bool = 1; 
plot_holo_rewards_bool = 1;
plot_neural_rewards_bool = 1;
plot_E1_bool = 1;
plot_E2_bool = 1;
plot_indirect_bool = 0; 

if(plot_x_in_time)
    plot_x = x_data*frame2time; 
else
    plot_x = x_data;
end

E1_data = neural(E1_chan, :);
E2_data = neural(E2_chan, :);
indirect_data = neural(indirect_chan, :);

E1_color = {};
for i=1:num_E1
    E1_color{i} = bmi_colors.E1;
end
E2_color = {};
for i=1:num_E2
    E2_color{i} = bmi_colors.E2;
end
indirect_color = {};
for i=1:num_indirect
    indirect_color{i} = bmi_colors.indirect;
end

if(plot_E1_bool && plot_E2_bool && plot_indirect_bool)
    traces = [indirect_data; E2_data; E1_data]; 
    trace_color = {indirect_color{:} E2_color{:} E1_color{:}};    
elseif(plot_E1_bool && plot_E2_bool)
    traces = [E2_data; E1_data]; 
    trace_color = {E2_color{:} E1_color{:}};
elseif(plot_E1_bool)
    traces = E1_data;
    trace_color = E1_color; 
elseif(plot_E2_bool)
    traces = E2_data; 
    trace_color = E2_color; 
elseif(plot_indirect_bool)
    traces = indirect_data; 
    trace_color = indirect_color;    
end

%-------------------------------------------------------------------------
h = figure;
hold on

err_bool = 0; 
traces_err = []; 
[min_val, max_val] = ...
    plot_traces_on_same_fig(plot_x, traces, err_bool, traces_err, trace_color);

if(plot_holo_stim_bool)
    %Mark the holo stimulations: 
    for i=1:length(datainfo.HoloStim_frames)
        if(plot_x_in_time)
            plot([datainfo.HoloStim_frames(i) datainfo.HoloStim_frames(i)]*frame2time, [0 max_val], 'b');
        else
            plot([datainfo.HoloStim_frames(i) datainfo.HoloStim_frames(i)]*frame2time, [0 max_val], 'b');
        end
    end
end

if(plot_holo_rewards_bool)
    %Mark the vta stim post holo stim rewards: 
    for i=1:length(reward.VTAholo)
        if(plot_x_in_time)
            plot([reward.VTAholo(i) reward.VTAholo(i)]*frame2time, [0 max_val], 'y');
        else        
            plot([reward.VTAholo(i) reward.VTAholo(i)], [0 max_val], 'y');
        end
    end
end

if(plot_neural_rewards_bool)
    %Mark the self rewards: 
    for i=1:length(reward.VTAreal)
        if(plot_x_in_time)
            plot([reward.VTAreal(i) reward.VTAreal(i)]*frame2time, [0 max_val], 'g');
        else
            plot([reward.VTAreal(i) reward.VTAreal(i)], [0 max_val], 'g');
        end
    end
end



%%
%Analyze the distribution over time separating holo_stim and vta_stim
%Confirm with Nuria that: rewardFrames / reward.VTA contain the frames
%during which VTA stim was delivered.  

%Loop over VTAholo frames, find the closest preceding holo stim, 
stim_delay = -ones(length(reward.VTAholo),1);
for i = 1:length(reward.VTAholo)
    reward_i = reward.VTAholo(i);
    time_diff = (reward_i - datainfo.HoloStim_frames);
    time_diff_before = time_diff(time_diff > 0);
    closest_holo_i = min(time_diff_before);
    stim_delay(i) = closest_holo_i;
end

%%
%Units of frames:
edges = 0:1:20; 
N = histc(stim_delay, edges);
h = figure;
bar(edges,N,'histc')
xlabel('num frames delay'); 
title('Delay between holostim and vtastim (in frames)'); 

edges = 0:1:500; 
N = histc(stim_delay, edges);
h = figure;
bar(edges,N,'histc')
xlabel('num frames delay'); 
title('Delay between holostim and vtastim (in frames)'); 
%%
sec_per_frame = 1/datainfo.FPS;
%Units of frames:
edges = 0:1:20; 
N = histc(stim_delay, edges);

edges_time = edges*sec_per_frame;
h = figure;
bar(edges_time,N,'histc')
xlabel('num sec delay'); 
title('Delay between holostim and vtastim (in sec)'); 

edges = 0:1:500; 
edges_time = edges*sec_per_frame;
N = histc(stim_delay, edges);
h = figure;
bar(edges_time,N,'histc')
xlabel('num sec delay'); 
title('Delay between holostim and vtastim (in sec)'); 

%%
%PSTH locked to holostim
%Assign neural source: 
neural = neurons.Fraw;

center_times = Holo_frames;
win = [-100 100]; 
trial_data = time_series2data_mat(neural, center_times, win);
%
reshape_size    = [size(trial_data,1)*size(trial_data,2) size(trial_data, 3)];
result_size     = [size(trial_data,1) size(trial_data,2)];
[sem_result, mean_result] = ...
    sem(reshape(trial_data, reshape_size), ones(reshape_size)); 
psth_data = reshape(mean_result, result_size); 
psth_sem = reshape(sem_result, result_size); 

%--------------------------------------------------------------------------
E1_psth = psth_data(E1_chan,:); 
E1_sem  = psth_sem(E1_chan,:);
E1_color = {};
for i=1:size(E1_psth,1)
    E1_color{i} = 'r';
end

E2_psth = psth_data(E2_chan,:);
E2_sem  = psth_sem(E2_chan,:); 
E2_color = {};
for i=1:size(E2_psth,1)
    E2_color{i} = 'k';
end

BMI_psth        = [E2_psth; E1_psth]; 
BMI_psth_sem    = [E2_sem; E1_sem]; 
BMI_color       = {E2_color{:} E1_color{:}};

%
x_data = win(1):win(2);
traces = BMI_psth;
traces_err = BMI_psth_sem;
trace_color = BMI_color;

h = figure; 
hold on; 
[min_val, max_val] = ...
    plot_traces_on_same_fig(x_data, traces, err_bool, traces_err, trace_color)
plot([0 0], [min_val max_val], 'b')
title('Red: E1, Black: E2, Blue: time of holo stim'); 
xlabel('Frame (15.2 FPS)'); 

%%
%PSTH locked to VTAstim following holostim
%Assign neural source: 
neural = neurons.Fraw;

center_times = reward.VTAholo;
win = [-100 100]; 
trial_data = time_series2data_mat(neural, center_times, win);
%
reshape_size    = [size(trial_data,1)*size(trial_data,2) size(trial_data, 3)];
result_size     = [size(trial_data,1) size(trial_data,2)];
[sem_result, mean_result] = ...
    sem(reshape(trial_data, reshape_size), ones(reshape_size)); 
psth_data = reshape(mean_result, result_size); 
psth_sem = reshape(sem_result, result_size); 

%--------------------------------------------------------------------------
E1_psth = psth_data(E1_chan,:); 
E1_sem  = psth_sem(E1_chan,:);
E1_color = {};
for i=1:size(E1_psth,1)
    E1_color{i} = 'r';
end

E2_psth = psth_data(E2_chan,:);
E2_sem  = psth_sem(E2_chan,:); 
E2_color = {};
for i=1:size(E2_psth,1)
    E2_color{i} = 'k';
end

BMI_psth        = [E1_psth; E2_psth]; 
BMI_psth_sem    = [E1_sem; E2_sem]; 
BMI_color       = {E1_color{:} E2_color{:}};

%
x_data = win(1):win(2);
traces = BMI_psth;
traces_err = BMI_psth_sem;
trace_color = BMI_color;

h = figure; 
hold on; 
[min_val, max_val] = ...
    plot_traces_on_same_fig(x_data, traces, err_bool, traces_err, trace_color)
plot([0 0], [min_val max_val], 'y')
title('Red: E1, Black: E2, Yellow: time of VTA stim following a holostim'); 
xlabel('Frame (15.2 FPS)'); 

%%
%PSTH locked to Self-Achieved VTA
center_times = reward.VTAreal;
win = [-100 100]; 
trial_data = time_series2data_mat(neural, center_times, win);

reshape_size    = [size(trial_data,1)*size(trial_data,2) size(trial_data, 3)];
result_size     = [size(trial_data,1) size(trial_data,2)];
[sem_result, mean_result] = ...
    sem(reshape(trial_data, reshape_size), ones(reshape_size)); 
psth_data = reshape(mean_result, result_size); 
psth_sem = reshape(sem_result, result_size); 
%--------------------------------------------------------------------------
E1_psth = psth_data(E1_chan,:); 
E1_sem  = psth_sem(E1_chan,:);
E1_color = {};
for i=1:size(E1_psth,1)
    E1_color{i} = bmi_colors.E1;
end

E2_psth = psth_data(E2_chan,:);
E2_sem  = psth_sem(E2_chan,:); 
E2_color = {};
for i=1:size(E2_psth,1)
    E2_color{i} = bmi_colors.E2;
end

BMI_psth        = [E1_psth; E2_psth]; 
BMI_psth_sem    = [E1_sem; E2_sem]; 
BMI_color       = {E1_color{:} E2_color{:}};

%
x_data = win(1):win(2);
traces = BMI_psth;
err_bool = 1; 
traces_err = BMI_psth_sem;
trace_color = BMI_color;

h = figure; 
hold on;
[min_val, max_val] = ...
    plot_traces_on_same_fig(x_data, traces, err_bool, traces_err, trace_color)
plot([0 0], [min_val max_val], 'g')
title('Red: E1, Black: E2, Green: time of self-achieved VTA stim'); 
xlabel('Frame (15.2 FPS)'); 

%%
%Do Factor Analysis on baseline 
%   BL_frames: [1x15567 double]
% Post_frames: [1x15266 double]

pre_data = neural([E1_chan; E2_chan], BL_frames); 
size(pre_data)

post_data = neural([E1_chan; E2_chan], Post_frames); 
size(post_data)

%%
%cov: col are var, row are obs
% C_pre = cov(pre_data.'); 
% C_post = cov(post_data.'); 
C_pre = corrcoef(pre_data.'); 
C_post = corrcoef(post_data.'); 

h = figure;
imagesc(C_pre); 
caxis([-.05 0.5]); 
colorbar
axis square
title('pre VTA stim training, stim neurons 1-4'); 

h = figure;
imagesc(C_post);
caxis([-.05 0.5]); 
colorbar
axis square
title('post VTA stim training, stim neurons 1-4'); 

%%
zDim = 3; 
X = pre_data*10000; 
[estParams, LL] = fastfa(X, zDim);
%ASSIGN:
estParams_pre = estParams; 
LL_pre = LL; 

% %
% zDim = 3; 
% X = post_data; 
% [estParams, LL] = fastfa(X, zDim);
% %ASSIGN:
% estParams_post = estParams; 
% LL_post = LL; 
% %%

%%
m = zDim; 
[lambda,psi,T,stats,F] = factoran(X,m)


