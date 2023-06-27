function c=intr(f,k,b1,b2)
t=cheb(b1,b2,k);
v=zeros(1,k);
for i=1:k
v(i)=f(t(i));
end
a=chebint(v,b1,b2,k);
c=intsr(a);
c=c*(b2-b1)/2;
end