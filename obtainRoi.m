function unitVals=obtainRoi(Im, neuronMask, com, units, sizNeur)
%{
Function to obtain the activity of each neuron, given a spatial filter
units -> unit which units to return
Im => Image 
neuronMask -> matrix for spatial filters with px*py*unit
com => position of the neurons
units => index of the neurons in the neuronMask
sizNeur => it has to be bigger than the radio of the neuron
%}

    if nargin < 5
        sizNeur = 20;
    end
    if nargin < 4
        units = 1:size(neuronMask,3);
    end

    unitVals = zeros(length(units),1);
    
    for auxu = 1: length(units)
        u = units(auxu);
        posmaxx = min([(com(u,1)+sizNeur), size(neuronMask,1)]);
        posminx = max([(com(u,1)-sizNeur), 1]);
        posmaxy = min([(com(u,2)+sizNeur), size(neuronMask,2)]);
        posminy = max([(com(u,2)-sizNeur), 1]);
        Imd = double(Im(posminy:posmaxy,posminx:posmaxx));
        auxmask = neuronMask(posminy:posmaxy,posminx:posmaxx,u);
        unitVals(auxu) = nansum(nansum(Imd.* auxmask/sum(auxmask)));
    end
end