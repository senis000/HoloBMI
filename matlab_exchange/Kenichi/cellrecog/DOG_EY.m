function out = DOG_EY(r)
r=r*2;
ker1 = DOG_N(r);
ker2 = DOG_N(r/2);

a=zeros(size(ker1));

s1=size(ker1,1);
s2=size(ker2,1);

a(round(s2/2):s1-round(s2/2)+1,round(s2/2):s1-round(s2/2)+1)=ker2;

out=ker1-a;