function [holoMask, Im] = makeMasksPrairie(channel)
    %{ Function to obtain masks in the red channel%}

    %channel = 1 ; %for the red channel
    %% connection to Prairie
    pl = actxserver('PrairieLink.Application');
    pl.Connect();

    % Prairie variables
    px = pl.PixelsPerLine();
    py = pl.LinesPerFrame();

    %take image from channel 1
    Im = pl.GetImage_2(channel, px, py);

    %% extract rois
    [mask, ~] = imFindCellsTM (Im, 11, 0.7, 9, 1, 0); %parameters depend on each image
    %displays the mask to see if we agree
    findCenter (mask, Im);
    holoMask = bwlabel(mask);
    pl.Disconnect();
end