
%%
mat     = zeros(512,512); 
temp    = zeros(512,512); 
dim     = 9; 
% iter    = 3; 
[mask, x_cen, y_cent] = addcell(mat,temp,dim);

%%
[x,y]=find_center_EY(mat,mask);

%%
[x,y]=findCenter(mask, mat)

%%
num_cells = max(max(mask)); 
for i =1:num_cells
    y(i)= round(mean(find(mask'==i))/size(mask,2));
    x(i)= round(mean(find(mask==i))/size(mask,1));
end