function [x,y]=findCenter(temp, mat, toplot)
%{
Function to find the center of mass for the detected cells
temp is the template
mat is the raw image
toplot is a flag to allow plotting
%}
    if nargin < 2
        toplot = false;
        mat = 0;
    elseif nargin < 3
        toplot = true;
    end

    y=zeros(max(max(temp)),1);
    x=zeros(max(max(temp)),1);
    for i=1:max(max(temp))
        y(i)= round(mean(find(temp'==i))/size(temp,2));
        x(i)= round(mean(find(temp==i))/size(temp,1));
    end

    if toplot
        figure ,
        if ~isempty(y)
            subplot (1,2,1)
            imagesc(temp), colormap bone, caxis([0 1]), hold on, scatter (x,fliplr(y), 'filled', 'r'), hold off

            subplot (1,2,2)
            imagesc(mat), colormap bone, caxis([-0 nanmean(nanmean(mat))*4]), hold on, scatter (x,fliplr(y), 'filled', 'r'), hold off
        else
            subplot (1,2,1)
            imagesc(temp), colormap bone, caxis([0 1])

            subplot (1,2,2)
            imagesc(mat), colormap bone, caxis([-0 nanmean(nanmean(mat))*4])

        end
    end
end