function out = imFilter2 (ker, in)

% ker: filter kernel
% in: input 2-D image
% our: output 2-D image
%
% This 2-dimensional filter routine is a modification of filter2.
% Around edges, filter2 calculates values assuming all values outside image are zero.
% So, if the pixel values does not appoach zero around the edge, 
% filter2 will underestimate the value.
%
% Kenichi Ohki 08/11/04

[xDim, yDim] = size(in);
temp = ones(xDim, yDim);
out = filter2(ker, in);
temp = filter2(ker, temp);
out = out./temp;