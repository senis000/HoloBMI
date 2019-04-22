function [flagBMI, flagVTAtrig, flagHolosched, flagVTAsched] = ...
    expt2bmi_flags(expt_str)
%Function which maps desired experiment type to experiment flags: 
%     0) BMI
%     flagBMI = true; (use self-generated hits to perform actions 
%         for now actions are just to possibly send VTA stim
%         in future, actions can include sending task feedback signal
%     flagHolo = false; flagVTAsched = false;
% 1) Pre-train (Holo -> VTA)
% flagVTAtrig = true; flagHolosched = true;
% This follows the vectorHolo schedule. If target pattern achieved,
% deliver VTA.
% flagBMI = false; flagVTAsched = false;
% 
% 2) No pre-training
% flagBMI = true; flagVTAtrig = false; flagHolosched=false;
% set experiment length to be length of pretrain + BMI
% flagVTAsched = false;
% 
% 3) Pre-train E3, Test E2
% baseline calibrate E3-E1, and E2-E1.  
% Pre-train E3: 
% flagVTAtrig = true; flagHolosched = true; flagBMI = false;
% flagVTAsched = false;
% 
% Test E2:
% flagVTAtrig = true; flagHolosched = false; flagBMI = true;
% flagVTAsched = false;
% 
% 4) Pre-train orthogonal
% baseline calibrate E2-E1 shuffle, and E2-E1
% Pretrain E2-E1 shuffle:
% flagVTAtrig = true; flagHolosched = true; flagBMI = false;
% flagVTAsched = false;
% 
% Test E2-E1: 
% flagVTAtrig = true; flagHolosched = false; flagBMI = true;
% flagVTAsched = false;
% 
% 5) Pre-train holo only
% flagVTA = false; flagHolo = true; flagBMI = false; 
% flagVTAsched = false;
% 
% 6) Random VTA
% flagVTAsched = true; flagVTAtrig = false; 
% flagBMI = false; flagHolosched = false;

expt_cell = {...
    'BMI', ...
    'HoloVTA_pretrain', ...
    'Holo_pretrain', ...
    'VTA_pretrain'}; 

if(strcmp(expt_str, 'BMI'))
    flagBMI         = true;
    flagVTAtrig     = true;
    flagHolosched   = false;
    flagVTAsched    = false;     
elseif(strcmp(expt_str, 'HoloVTA_pretrain'))
    flagBMI         = false;
    flagVTAtrig     = true;
    flagHolosched   = true;
    flagVTAsched    = false; 
elseif(strcmp(expt_str, 'Holo_pretrain'))
    flagBMI         = false;
    flagVTAtrig     = false;
    flagHolosched   = true;
    flagVTAsched    = false;    
elseif(strcmp(expt_str, 'VTA_pretrain'))
    flagBMI         = false;
    flagVTAtrig     = false;
    flagHolosched   = false;
    flagVTAsched    = true;      
else
    disp('You did not enter an acceptable expt str')
end