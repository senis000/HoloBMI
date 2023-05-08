function [mask] = create_mask_rod2 (imag,dim)

%[mask] = addcell (mat,mask,dim)
% imag => the original raw image
% mask => the mask to modify
% dim => diameter of the cell (even number)


center = floor(dim/2)+1;
radius = round (dim/2);

mask = zeros(size(imag));

while max(max(imag))> 100
    [~, new_cells] = CircularHough_Grd(double(imag), [2 15],1000,4);

    new_cells = round(new_cells);
    num_cell = max(max(mask)) + 1;
    
    for k = 1: size(new_cells,1)
        for i = 1:dim
            for j = 1:dim
               if (norm([i,j]-[center,center])<radius) && (new_cells(k,2)+round(i-dim/2)<= size(mask,1)) && (new_cells(k,1)+round(j-dim/2) <= size(mask,2)) ...
                       && (new_cells(k,2)+round(i-dim/2) >= 1) && (new_cells(k,1)+round(j-dim/2) >= 1)
                    mask(new_cells(k,2)+round(i-dim/2),new_cells(k,1)+round(j-dim/2)) = num_cell;
               end
            end
        end
        num_cell = num_cell + 1;
    end
    
    

    imag(imag>max(max(imag))/3)= 0 ;
end


[x_cen,y_cent]=find_center_EY(imag,mask);