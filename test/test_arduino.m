% test_arduino.m
%%
%Find way to do pure tones + reward sound + reward delivery with 

%%

a = arduino('COM12', 'Mega2560')
disp('arduino complete')
%%
tic
pin = 'D3'
frequency = 1000
duration = 0.3
playTone(a,pin,frequency,duration)
toc
%%
duration = 0.1
%%
% freq_sweep = 1000:500:16000;
freq_sweep = 6000:500:19000;

freq_sweep = 6000:500:19000;

% h = figure;
% plot(freq_sweep)
for j = 1
    for i=1:length(freq_sweep)
        tic
        freq_i = freq_sweep(i); 
        playTone(a,pin,freq_i,duration);
        toc
        pause(0.03)
    end
end

%%
freq_sweep = 1000:500:16000;
% h = figure;
% plot(freq_sweep)
for j = 1
    for i=1:length(freq_sweep)
        tic
        freq_i = freq_sweep(length(freq_sweep)-i+1); 
        playTone(a,pin,freq_i,duration);
        toc
        pause(0.03)
    end
end

%%


%%
duration = 1
playTone(a,pin,5000,duration);

%%
%Compare to matlab sound: 
Fs = 10000; 
xrnd = randn(Fs*1,1);
reward_sound = audioplayer(xrnd, Fs); %Play sound using: play()

%%
tic
sound(1e5*xrnd, Fs)
toc

%%
t=1:Fs;
f = 1000; 
y=cos(2*pi/f*t);
% h = figure;
% plot(y)
tic
soundsc(y, Fs)
toc

%%
tic
soundsc(xrnd, Fs)
toc


%%

%%
% writeDigitalPin(a, 'D13', 0);
% pause(1)
tic
writeDigitalPin(a, 'D13', 1);
pause(0.135)
writeDigitalPin(a, 'D13', 0);
toc
%Flipping on and off D13 costs 20-25ms...

%%
%Compare to turning on and off nidaq

%%
% Equation for figuring out the maximum number of 1/4 octave increments for
% a given range
y = 20000;  % final frequency
a = 1000;   % initial frequency
b = 1.25;   % increment

max_inc = (log(y)-log(a))/log(b);

% Now to Generate the desired frequencies
freqs   = zeros(floor(max_inc),1);

for n = 1:floor(max_inc)
freqs(n) = a*(b^n);
end

%%
% Set sampling rate
Fs = 44100;

% Set time duration (100ms per tone currently)
x = 0:Fs;

% Construct Gaussian Window
sig = 44100/(2*pi);

% Hamming window the signal to make transitions smooth
ham = 0.54 - 0.46.*cos(2*pi*(x/length(x)));

for n = 1:length(freqs)
    y=ham.*sin(2*pi*(freqs(n)/Fs).*x);
    plot(x,y)
    pause
    soundsc([y' y'],Fs) % sound([y' y'],44100)
end

%%
y=ham.*sin(2*pi*(freqs(5)/Fs).*x)
soundsc(y, Fs)

%%
freq_i = 7000; 
x = 0:Fs*10;
y=sin(2*pi*(freq_i/Fs).*x);
soundsc(y, Fs)

%%
duration = 1
playTone(a,pin,freq_i,duration);

%%
max_freq    = 19000;
min_freq    = 0; 
freq_range  = max_freq - min_freq; 

len = 200;
rndTone = freq_range*rand(len,1) + min_freq; 

for i=1:len
%     tic
    playTone(a, pin, rndTone(i), 0.0005); 
%     toc
end


