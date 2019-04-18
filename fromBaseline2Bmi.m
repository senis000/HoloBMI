function [E1BMI, E2BMI] = fromBaseline2Bmi(E1, E2)
%{
 Function that does the transition for variables from Baseline to BMI
E1,E2 => Index in the spatial filter matrix
%}

    savePath = ['F:/VivekNuria/', animal, '/',  day, '/'];
    if ~exist(savePath, 'dir')
            mkdir(savePath);
    end
    load([savePath, 'redcomp.mat'], 'AComp');
    ensemble = [E1, E2];
    AComp = AComp(:, ensemble);
    E1BMI = 1:length(E1);
    E2BMI = length(E1)+1:length(E1)+length(E2);
    save(savePath + "redcompBMI.mat", 'AComp')
        
end