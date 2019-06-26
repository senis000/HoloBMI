function strcMask=obtainStrcMaskfromMask(mask)
%{
Function to create a structure with a reduce image given a spatial filter 
AComp => Spatial filter
px,py => shapes of the image
strcMask -> structure with the matrix for spatial filters with px*py*unit
and the positions of that mask
units => index of the neurons in the neuronMask
%}

    for u = 1: max(max(mask))
        u
        auxMask = mask;
        auxMask(auxMask~=u) = 0;
        posx = find(sum(auxMask,1)~=0);
        posy = find(sum(auxMask,2)~=0);
        strcMask.maxx(u) = posx(end);
        strcMask.minx(u) = posx(1);
        strcMask.maxy(u) = posy(end);
        strcMask.miny(u) = posy(1);
        strcMask.neuronMask{u} = auxMask(posy(1):posy(end), posx(1):posx(end));
    end
end