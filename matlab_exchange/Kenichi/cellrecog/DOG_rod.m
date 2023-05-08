function ker = DOG_rod(r)

r=r*5.1;
r1=r/5;         % zero-crossing at r
r2=r1*2;

dim=ceil(r2*3);
ker=fspecial('log',dim,r1)-(fspecial('log',dim,r2)*0.6)*15;
return;