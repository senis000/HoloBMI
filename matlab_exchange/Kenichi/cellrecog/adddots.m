function [mask] = adddots (X,Y,mask,dim)


new_dot = [X;Y]';
center = floor(dim/2)+1;
radius = round (dim/2);

for i = 1:length(X)
    mask(Y(i),X(i)) = i;
end
num_cell = max(max(mask)) + 1;

for k = 1: size(new_dot,1)
    for i = 1:dim
        for j = 1:dim
           if (norm([i,j]-[center,center])<radius) && (new_dot(k,2)+round(i-dim/2)<= size(mask,1)) && (new_dot(k,1)+round(j-dim/2) <= size(mask,2))
                mask(new_dot(k,2)+round(i-dim/2),new_dot(k,1)+round(j-dim/2)) = num_cell;
           end
        end
    end
    num_cell = num_cell + 1;
end
