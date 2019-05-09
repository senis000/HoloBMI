function [holoMask, Im, Img, px, py] = makeMasksPrairie()
    %{ Function to obtain masks in the red channel%}

    %channel = 1 ; %for the red channel
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
    [mask, ~] = imFindCellsTM (Im-Img/10, 11, 0.5, 9, 1, 0); %parameters depend on each image
    %displays the mask to see if we agree
    findCenter (mask, Im);
    holoMask = bwlabel(mask);
    pl.Disconnect();
end