


clc
close all; 

im = zeros(512,512); 
h = figure;
imagesc(im); 
axis square

r = round(ginput(1))

hold on;
scatter(r(1),r(2))

%Column, Row
%X, 512-Y+1
%

% Column    = X
% Row       = 512-Y+1
% X         = Column
% Y         = 512-Row+1
%This means if you want to scatter (x,y) data onto an image
%

%Top Left: 1,1
%Top Right: 512, 1
%Bot Right: 512,512
%Bot Left: 1, 512

%%
drawRois(numArea);  %draw rois

