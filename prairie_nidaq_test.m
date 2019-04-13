% Connects prairie and nidaq to test for the speed of the closed loop

%% Variables
chan_idx = 2;
sleep_time = 0;
ttl_frame = 15;
buffer_size = 1048576;
number_frames = 100; 
gi_bool = 1 ;
rrds_bool = 0;


%% Prepare for Prairie
% connection to Prairie
pl = actxserver('PrairieLink.Application');
pl.Connect()

% Prairie variables
px = pl.PixelsPerLine();
py = pl.LinesPerFrame();
pixels_per_frame = px*py;
spp = pl.SamplesPerPixel();

% Prairie commands
pl.SendScriptCommands("-srd True 0");
pl.SendScriptCommands("-lbs True 0");


%% Prepare the nidaq
s = daq.createSession('ni');
addDigitalChannel(s,'dev5','Port0/Line0:0','OutputOnly');

%% Run the test

if gi_bool
    data_img = zeros(number_frames, px,py) + nan;
end
if rrds_bool
    data_rrds = zeros(number_frames, px,py) + nan;
    index_rrds = 1;
    samplesRead = pl.ReadRawDataStream(index_rrds);  % to clean the buffer
end

%start the live scan
pl.SendScriptCommands("-lv on");    

time_clock_gi = zeros(number_frames+1,1);
time_clock_rrds = zeros(number_frames+1,1);

tic;

time_clock_rrds(1) = toc;
time_clock_gi(1) = toc;


for i=1:number_frames
        %On this condition, trigger the LED on
    if i==ttl_frame
        outputSingleScan(s,1);
        disp('TTL OUT')
    else
        outputSingleScan(s,0);
    end

    if gi_bool
        Im = pl.GetImage_2(chan_idx, px, py);
        time_clock_gi(i+1) = toc;
    end
    if rrds_bool
            samplesRead = pl.ReadRawDataStream(px*py*spp);
        time_clock_rrds(i+1) = toc;
    end
        
    if gi_bool
        data_img(i,:,:) = Im;
    end
    if rrds_bool
        aux_samplesRead = samplesRead(1:px*py*spp);
        datares = reshape(aux_samplesRead, [px, py]);
        for j= 1:2:py
            datares(j,:) = datares(j,end:-1:1);
        end
        data_rrds(i, :, :) = datares; 
    end

    pause(sleep_time);
end
pl.SendScriptCommands("-lv off");  
disp('done!');
