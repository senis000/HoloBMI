function strcMask=obtainStrcMask(AComp, px, py)
%{
Function to obtain the activity of each neuron, given a spatial filter
units -> unit which units to return
Im => Image 
strcMask -> structure with the matrix for spatial filters with px*py*unit
and the positions of that mask
units => index of the neurons in the neuronMask
%}

    for u = 1: size(AComp,2)
        a1 = reshape(full(AComp(:,u)),px,py);
        posx = find(sum(a1,1)~=0);
        posy = find(sum(a1,2)~=0);
        strcMask.maxx(u) = posx(end);
        strcMask.minx(u) = posx(1);
        strcMask.maxy(u) = posy(end);
        strcMask.miny(u) = posy(1);
        strcMask.neuronMask{u} = a1(posy(1):posy(end), posx(1):posy(end));
    end
end