function a=chebint(f,b1,b2,k)
a=zeros(1,k);
t=cheb(-1,1,k);
for i=2:k
a(i)=2/k*sum(f.*cos((i-1)*acos(t)));
end
a(1)=1/k*sum(f);
end