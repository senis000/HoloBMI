% test_arduino.m

a = arduino('COM10', 'Mega2560')

%%
writeDigitalPin(a, 'D13', 0);
pause(1)
tic
writeDigitalPin(a, 'D13', 1);
pause(0.001)
writeDigitalPin(a, 'D13', 0);
toc
%Flipping on and off D13 costs 20-25ms...

%%
%Compare to turning on and off nidaq