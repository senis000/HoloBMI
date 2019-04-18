function strcMask=obtainStrcMask(AComp, px, py)
%{
Function to create a structure with a reduce image given a spatial filter 
AComp => Spatial filter
px,py => shapes of the image
strcMask -> structure with the matrix for spatial filters with px*py*unit
and the positions of that mask
units => index of the neurons in the neuronMask
%}

    for u = 1: size(AComp,2)
        a1 = reshape(full(AComp(:,u)),py,px);
        posx = find(sum(a1,1)~=0);
        posy = find(sum(a1,2)~=0);
        strcMask.maxx(u) = posx(end);
        strcMask.minx(u) = posx(1);
        strcMask.maxy(u) = posy(end);
        strcMask.miny(u) = posy(1);
        strcMask.neuronMask{u} = a1(posy(1):posy(end), posx(1):posx(end));
    end
end