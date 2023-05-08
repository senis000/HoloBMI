function [stack,info] = readstack(id, Zs, Ts, Cs, nbinning, period, runs)

% open a nikon nd2 or a zeiss lsm file
% 
% id: filename, or path\filename
%
% Zs: Z to read
% Ts: T to read
% Cs: C to read
%
% if Zs, Ts, Cs are not specified or empty, read all frames
%
% nbinning: bin size (XY), make images smaller
%   if empty or 1, no binning
%
% period: if you want to average data over repeats, 
%   specify period of one repeat; if empty, no average
% runs: which runs you want to average. if empty, use all repeats.
%
% e.g. [stack, info] = =readstack3(fname,[],[101:200],1); 
% read only 1st channel, between T=101~200
%
% stack: output stack, up to 5 dimension (XYZTC)
%
% info: information of stack
%
% modified from bfopen.m to read stack data in 5 dimension matrix (XYZTC).
%
%   2011. 5. 11. Kenichi Ohki
%
% A script for opening microscopy images in MATLAB using Bio-Formats.
%
% The function returns a list of image series; i.e., a cell array of cell
% arrays of (matrix, label) pairs, with each matrix representing a single
% image plane, and each inner list of matrices representing an image
% series. See below for examples of usage.
%
% Portions of this code were adapted from:
% http://www.mathworks.com/support/solutions/en/data/1-2WPAYR/
%
% This method is ~1.5x-2.5x slower than Bio-Formats's command line
% showinf tool (MATLAB 7.0.4.365 R14 SP2 vs. java 1.6.0_20),
% due to overhead from copying arrays.
%
% Thanks to all who offered suggestions and improvements:
%     * Ville Rantanen
%     * Brett Shoelson
%     * Martin Offterdinger
%     * Tony Collins
%     * Cris Luengo
%     * Arnon Lieber
%     * Jimmy Fong
%
% NB: Internet Explorer sometimes erroneously renames the Bio-Formats library
%     to loci_tools.zip. If this happens, rename it back to loci_tools.jar.
%
%



tic

%javaaddpath('D:/public/matlabcodes/ohkilab_official_programs/java/loci_tools.jar');

r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);

r.setId(id);
r.setSeries(0);

info = getBFfileInfo(r);

if nargin < 2 || isempty(Zs)
    Zs = [1:info.Z];
end
if nargin < 3 || isempty(Ts)
    Ts = [1:info.T];
end
if nargin < 4 || isempty(Cs)
    Cs = [1:info.C];
end

if nargin < 5 || isempty(nbinning)
    nbinning =1;
end

if nargin < 6
    period =[];
end

if nargin < 7 && ~isempty(period)
    runs = [1:floor(length(Ts)/period)];
end

if ~isempty(find(Zs<1)) || ~isempty(find(Zs>info.Z)) 
    display('Error: Z is out of range');
    stack=[];
    return;
end

if ~isempty(find(Ts<1)) || ~isempty(find(Ts>info.T)) 
    display('Error: T is out of range');
    stack=[];
    return;
end

if ~isempty(find(Cs<1)) || ~isempty(find(Cs>info.C)) 
    display('Error: C is out of range');
    stack=[];
    return;
end

Xsize=floor(info.X/nbinning);
Ysize=floor(info.Y/nbinning);

if isempty(period)
    
    stack=zeros(Xsize, Ysize, length(Zs), length(Ts), length(Cs), info.type);

    fprintf('Reading \n    ');

    for t = 1:length(Ts)
        for z = 1:length(Zs)
            for c = 1:length(Cs)
                i = Cs(c) + (Zs(z)-1)*info.C + (Ts(t)-1)*info.C*info.Z;
                n = c+(z-1)*length(Cs)+(t-1)*length(Cs)*length(Zs);
                if mod(n, 7200) == 0
                    fprintf('\n    ');
                end
                if mod(n, 100)==0
                    fprintf('.');
                end
                if nbinning == 1
                    stack(:,:,z,t,c) = readBFimage(r, i);
                else
                    stack(:,:,z,t,c) = binning(readBFimage(r, i), [nbinning nbinning]);
                end
            end
        end
    end
else
    Nrep=floor(length(Ts)/period);
    if ~isempty(find(runs<1)) || ~isempty(find(runs>Nrep)) || Nrep<1
        display('Error: runs are out of range');
        stack=[];
        return;
    end    
        
    stack=zeros(Xsize, Ysize, length(Zs), period, length(Cs));

    fprintf('Reading \n    ');
    for run = runs
        for t = 1:period
            for z = 1:length(Zs)
                for c = 1:length(Cs)
                    i = Cs(c) + (Zs(z)-1)*info.C + (Ts(t+(run-1)*period)-1)*info.C*info.Z;
                    n = c+(z-1)*length(Cs)+(t+(run-1)*period-1)*length(Cs)*length(Zs);
                    if mod(n, 7200) == 0
                        fprintf('\n    ');
                    end
                    if mod(n, 100)==0
                        fprintf('.');
                    end
                    if nbinning == 1
                        stack(:,:,z,t,c) = stack(:,:,z,t,c)+double(readBFimage(r, i));
                    else
                        stack(:,:,z,t,c) = stack(:,:,z,t,c)+double(binning(readBFimage(r, i), [nbinning nbinning]));
                    end
                end
            end
        end
    end
end

fprintf('\n');

r.close();

toc
end

