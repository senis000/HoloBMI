% Main protocol for the experiment

% Select area of interest with the 2p

% find red neurons
[mask, Im, px, py] = makeMasksPrairie(channel);

% delete neurons that we don't want by position on image
mask = deleteMask(Im, mask, numArea);

% HOLO STIM with the mask of the red components as input!!!!

% obtain A from onacid and bring it to matlab (check python file
% "obtain_componenents" -> 
%         'redLabel' : label as True for components marked as red,
%         'indRed' : index of such neurons,
%         'AComp' : sparse matrix with the spatial filters
fname = [folder , 'redcom.mat'];

load(fname)
numberNeurons = size(AComp,2);
AFull = reshape(full(AComp), px, py, numberNeurons);
%TODO do we want to remove neurons here? 

% acquire baseline 