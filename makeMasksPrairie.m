function [holoMask, Im, Img, px, py] = makeMasksPrairie(thres)
    %{ 
Function to obtain masks in the red channel
thres -> 0-1 high - strict threshold for the pattern recognition
%}
 
    %channel = 1 ; %for the red channel
    if nargin < 1
        thres = 0.5;
    end
        
    
    %% connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect();
    disp('Connecting to prairie')
    pause(2);

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    %take image from channel 1
    Im = pl.GetImage_2(1, px, py);
    Img = pl.GetImage_2(2, px, py);

    %% extract rois
    [mask, ~] = imFindCellsTM (Im-Img/10, 9, thres, 7, 1, 0); %parameters depend on each image
    %displays the mask to see if we agree
    findCenter (mask, Im);
    holoMask = bwlabel(mask);
    pl.Disconnect();
end