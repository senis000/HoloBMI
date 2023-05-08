function [stack, param] = readlsm(fname)

% open a zeiss LSM file
% using Francois Nedelec's tiffread.m
%
% fname: filename, or path\filename
%
% stack: output stack, up to 5 dimension (XYZTC)
%
% param: parameters of the file
%
%   2011. 4. 29.    Kenichi Ohki

first=tiffread(fname,1);

param=first.lsm;

Nchannel=param.DimensionChannels;

Nx=param.DimensionX;
Ny=param.DimensionY;
Nz=param.DimensionZ;
Nt=param.DimensionTime;

data=tiffread(fname,1:2:Nz*Nt*2-1);
% skip every other frames. Thumbnails are saved in even frames.

stack=zeros(Nx,Ny,Nz*Nt,Nchannel,class(first.data{1}));

for channel=1:Nchannel
    for i=1:Nz*Nt
       stack(:,:,i,channel)=data(i).data{channel}';
    end
end

stack=reshape(stack,Nx,Ny,Nz,Nt,Nchannel);




