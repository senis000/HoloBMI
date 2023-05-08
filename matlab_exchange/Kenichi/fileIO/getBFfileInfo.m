
function info = getBFfileInfo(r)

% r: java object of Bio-Formats, or file name
%
% extracted from bfopen.m  by Kenichi Ohki 2011.5.16.

if ~isjava(r)
    if ischar(r)
        fname=r;
        clear r
        r = loci.formats.ChannelFiller();
        r = loci.formats.ChannelSeparator(r);
        r.setId(fname);
        r.setSeries(0);
    else
        display('input should be java object or file name');
        info=[];
        return
    end
end

info.X = r.getSizeX();
info.Y = r.getSizeY();
info.Z = r.getSizeZ();
info.T = r.getSizeT();
info.C = r.getSizeC();

pixelType = r.getPixelType();
info.pixelType = pixelType;
info.bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
info.fp = loci.formats.FormatTools.isFloatingPoint(pixelType);
info.sgn = loci.formats.FormatTools.isSigned(pixelType);
info.bppMax = power(2, info.bpp * 8);
info.type = char(loci.formats.FormatTools.getPixelTypeString(pixelType));

info.little = r.isLittleEndian();

info.numImages = r.getImageCount();

info.Metadata = r.getMetadata();
end
