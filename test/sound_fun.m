% chunks of code for sound generation and fun auditory facts

% default sampling rate 8192 hz for "sound" command
% Dell  A525 speaker response bandwidth 45-20000 hz
% rat auditory range 1000-64000 hz
% rat auditory discrimination ~ 1/4 octave
% octave = frequency doubling

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

h = figure;
plot(freqs)

%%
% Set sampling rate
Fs = 44100;

% Set time duration (100ms per tone currently)
x = 0:Fs/10;

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