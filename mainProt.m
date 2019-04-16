% Main protocol for the experiment

% Select area of interest with the 2p

% find neurons
[mask, Im] = makeMasksPrairie(channel);

% delete neurons that we don't want by position on image
mask = deleteMask(Im, mask, numArea);

% obtain A from onacid and bring it to matlab
fname = [folder , 'redcom.mat'];

load(fname)
% numberNeurons = nanmax(nanmax(mask));
% neuronMask = zeros(numberNeurons, px, py);
% for i = 1:numberNeurons
%     neuronMask(i,:,:) = bwl==i;
% end

%TODO do we want to remove neurons here? 