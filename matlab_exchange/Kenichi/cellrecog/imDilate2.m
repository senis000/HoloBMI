function out=imDilate2(in, diam)

% imdilate2: dilate image with a circle of "diam"
%
% Kenichi Ohki  04.28.2008

radius=diam/2;
dim=ceil(radius)*2+1;
center=ceil(radius)+1;
NHOOD=zeros(dim);
for x=1:dim
    for y=1:dim
        if norm([x,y]-[center,center])<=radius
            NHOOD(x,y)=1;
        end
    end
end
out=imdilate(in,NHOOD);
    