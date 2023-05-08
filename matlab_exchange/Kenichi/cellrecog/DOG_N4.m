function ker = DOG_N4(r)
r=r*2.2;
r1=r/1.9227;         % zero-crossing at r
r2=r1*2;

dim=ceil(r2*2.5);
ker=fspecial('log',dim,r1)-(fspecial('log',dim,r2)*1.1)*15;
return;