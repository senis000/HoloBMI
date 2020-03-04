function [mask] =  circle_mask(dim, ctr_row, ctr_col, r)

mask = zeros(dim, dim); 
for row = 1:dim
    for col = 1:dim
        dist = norm([col row] - [ctr_col ctr_row]); 
        if dist <= r
            mask(row, col) = 1; 
        end
    end
end