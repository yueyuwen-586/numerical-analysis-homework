function [u,ud]=directode2solver(a,c,K,p0,q0,f0,A,gamma1,gamma2,nodelist)
%这是direct方法线性方程求解器，输入区间端点a、c，插值次数K，
%二阶方程参数函数p0、q0、f0，边值条件矩阵A=[zetal0,zetar0;zetal1,zetar1]，
%边界条件值gamma1,gamma2,分段插值的子区间端点序列nodelist。输出Func类函数u,ud作为u和u'的解。
%在directode2solvertest版本中额外输入：真实解solution；输出：误差。

l=max(size(nodelist))-1;
p=p0;
q=@(t) q0(t)-Q0(A);
b=decompose(A,a,c,gamma1,gamma2);
ui=@(t) b(2)+b(1).*t;
uid=@(t) b(1);
uidd=@(t) 0;
f=@(t) f0(t)-(uidd(t)+p0(t).*uid(t)+q0(t).*ui(t));
s=S(A,a,c);
g1d=@(t) gld(t,A,a,c);
g2d=@(t) grd(t,A,a,c);
g1=@(t) gl(t,A,a,c);
g2=@(t) gr(t,A,a,c);
psi1=@(t) psil(t,p,q,A,a,c);
psi2=@(t) psir(t,p,q,A,a,c);
sigma=globalPinverse(f,K,nodelist,g1,g2,psi1,psi2);
sigmafunc=Func(nodelist,sigma);
J1=zeros(1,l);
J2=zeros(1,l);
coefu=zeros(K,l);
coefud=zeros(K,l);
for i=1:l-1
    prodl=@(t) g1(t).*compute(sigmafunc,t);
    prodr=@(t) g2(t).*compute(sigmafunc,t);
    J1(i+1)=J1(i)+integral(prodl,nodelist(i),nodelist(i+1));
    J2(l-i)=J2(l-i+1)+integral(prodr,nodelist(l+1-i),nodelist(l+2-i));
end
for i=1:l
    [w,wd]=localsolveequationwithoutnode(nodelist(i),nodelist(i+1),sigma(:,i)',K,J1(i),J2(i),ui,uid,s,g1,g2,g1d,g2d);
    coefu(:,i)=w';
    coefud(:,i)=wd';
end
u=Func(nodelist,coefu);
ud=Func(nodelist,coefud);
end