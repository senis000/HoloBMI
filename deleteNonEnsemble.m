function EnsembleMask=deleteNonEnsemble (AComp, E2indx, px, py)
%{
Function to obtain the mask of only the ensemble neurons
AComp => Spatial filter
px,py => shapes of the image
E2indx => index of the neurons to keep in the mask 
%}

AEns = permute(reshape(full(AComp(:,E2indx)),py,px, length(E2indx)), [2,1,3]);
AEns(AEns>0.1) = 1;

EnsembleMask = sum(AEns,3);

se = strel('diamond',10);

EnsembleMask = bwlabel(imerode(imdilate(EnsembleMask, se), se));



end