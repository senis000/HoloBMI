function [mask, x_cen, y_cent] = addcell (mat,temp,dim,iter)

%[mask] = addcell (mat,temp,dim)
% mat => the original raw image
% temp => the mask to modify
% dim => diameter of the cell (even number)

if nargin <4
    iter = 0;
end

find_center_EY(mat,temp);
pause
new_cells = round(ginput);

%to add in more than one region, pause to focus in the area with the proper
%caxis
for i = 1:iter-1
    pause
    new_cells = [new_cells ; round(ginput)];
end

center = floor(dim/2)+1;
radius = round (dim/2);
num_cell = max(max(temp)) + 1;

for k = 1: size(new_cells,1)
    for i = 1:dim
        for j = 1:dim
           if (norm([i,j]-[center,center])<radius) && (new_cells(k,2)+round(i-dim/2)<= size(temp,1)) && (new_cells(k,1)+round(j-dim/2) <= size(temp,2))...
                    && (new_cells(k,2)+round(i-dim/2) >= 1) && (new_cells(k,1)+round(j-dim/2) >= 1)
                temp(new_cells(k,2)+round(i-dim/2),new_cells(k,1)+round(j-dim/2)) = num_cell;
           end
        end
    end
    num_cell = num_cell + 1;
end
mask = temp;

[x_cen,y_cent]=find_center_EY(mat,mask);