function a=Pinverse(f,K,b1,b2,g1,g2,psi1,psi2)
t=cheb(b1,b2,K);
A=zeros(K,K);
for j=1:K
  h1=@(t) g1(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
  h2=@(t) g2(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
  c1=intl(h1,K,b1,b2);
  c2=intr(h2,K,b1,b2);
  for i=1:K
    A(i,j)=cos((j-1)*acos(2*(t(i)-b1)/(b2-b1)-1))+psi1(t(i))*chebvalue(c1,b1,b2,t(i))+psi2(t(i))*chebvalue(c2,b1,b2,t(i));
  end
end
c=f(t);
a=A\(c');
a=a';
end