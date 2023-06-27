function [w,wd]=localsolveequationwithoutnode(b1,b2,sig,K,J1,J2,ui,uid,s,g1,g2,g1d,g2d)
u=zeros(1,K);
ud=zeros(1,K);
t=cheb(b1,b2,K);
h1=@(t) chebvalue(sig,b1,b2,t)*g1(t);
h2=@(t) chebvalue(sig,b1,b2,t)*g2(t);
c1=intl(h1,K,b1,b2);
c2=intr(h2,K,b1,b2);
for j=1:K
    u(j)=ui(t(j))+g2(t(j))*(J1+chebvalue(c1,b1,b2,t(j)))/s+g1(t(j))*(J2+chebvalue(c2,b1,b2,t(j)))/s;
    ud(j)=uid(t(j))+g2d(t(j))*(J1+chebvalue(c1,b1,b2,t(j)))/s+g1d(t(j))*(J2+chebvalue(c2,b1,b2,t(j)))/s;
end
w=chebint(u,b1,b2,K);
wd=chebint(ud,b1,b2,K);
end