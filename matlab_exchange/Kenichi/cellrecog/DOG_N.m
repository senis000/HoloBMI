function ker = DOG_N(r)
r1=r/1.9227;         % zero-crossing at r
r2=r1*2;

dim=ceil(r2*5);
ker=fspecial('gaussian',dim,r1)-fspecial('gaussian',dim,r2);
return;