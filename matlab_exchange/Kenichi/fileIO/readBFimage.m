
function arr = readBFimage(r, i)

% r: java object of Bio-Formats
% i: i-th frame is read. i starts from 1.
%
% extracted from bfopen.m  by Kenichi Ohki 2011.5.16.

% check Bio-Formats version, since makeDataArray2D function requires trunk
bioFormatsVersion = char(loci.formats.FormatTools.VERSION);
isBioFormatsTrunk = versionCheck(bioFormatsVersion, 5, 0);

canTypecast = versionCheck(version, 7, 1);

info = getBFfileInfo(r);

plane = r.openBytes(i - 1);
% convert byte array to MATLAB image
if isBioFormatsTrunk && (info.sgn || ~canTypecast)
    % can get the data directly to a matrix
    arr = loci.common.DataTools.makeDataArray2D(plane, ...
        info.bpp, info.fp, info.little, info.Y);
else
    % get the data as a vector, either because makeDataArray2D
    % is not available, or we need a vector for typecast
    arr = loci.common.DataTools.makeDataArray(plane, ...
        info.bpp, info.fp, info.little);
end

arr=typecastJava(arr,info.sgn,info.bppMax);
arr=reshape(arr, info.X, info.Y);

end       
         
function result = typecastJava (data, sgn, bppMax)
% Java does not have explicitly unsigned data types;
% hence, we must inform MATLAB when the data is unsigned

% check MATLAB version, since typecast function requires MATLAB 7.1+
canTypecast = versionCheck(version, 7, 1);

if ~sgn
    if canTypecast
        % TYPECAST requires at least MATLAB 7.1
        % NB: arr will always be a vector here
        switch class(data)
            case 'int8'
                result = typecast(data, 'uint8');
            case 'int16'
                result = typecast(data, 'uint16');
            case 'int32'
                result = typecast(data, 'uint32');
            case 'int64'
                result = typecast(data, 'uint64');
        end
    else
        % adjust apparent negative values to actual positive ones
        % NB: arr might be either a vector or a matrix here
        mask = data < 0;
        adjusted = data(mask) + bppMax / 2;
        switch class(data)
            case 'int8'
                result = uint8(data);
                adjusted = uint8(adjusted);
            case 'int16'
                result = uint16(data);
                adjusted = uint16(adjusted);
            case 'int32'
                result = uint32(data);
                adjusted = uint32(adjusted);
            case 'int64'
                result = uint64(data);
                adjusted = uint64(adjusted);
        end
        adjusted = adjusted + bppMax / 2;
        data(mask) = adjusted;
    end
end
end


