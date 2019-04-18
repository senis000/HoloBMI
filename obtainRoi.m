function unitVals=obtainRoi(Im, strcMask, units)
%{
Function to obtain the activity of each neuron, given a spatial filter
units -> unit which units to return
Im => Image 
strcMask -> structure with the matrix for spatial filters with px*py*unit
and the positions of that mask
units => index of the neurons in the neuronMask
%}

    if nargin < 3
        units = 1:size(strcMask.neuronMask,3);
    end

    unitVals = zeros(length(units),1);

    for auxu = 1: length(units)
        u = units(auxu);
        posmaxx = strcMask.maxx(u);
        posminx = strcMask.minx(u);
        posmaxy = strcMask.maxy(u);
        posminy = strcMask.mainy(u);
        Imd = double(Im(posminy:posmaxy,posminx:posmaxx));
        unitVals(auxu) = nansum(nansum(Imd.* strcMask.neuronMask{u}/nansum(strcMask.neuronMask{u})));
    end
end