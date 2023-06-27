function [u,ud]=ode2solver(a,c,K,p0,q0,f0,A,gamma1,gamma2,C,Tol)
%这是mesh refinement算法线性方程求解器，输入区间端点a、c，插值次数K，
%二阶方程参数函数p0、q0、f0，边值条件矩阵A=[zetal0,zetar0;zetal1,zetar1]，
%边界条件值gamma1,gamma2,算法中的C和Tol。输出Func类函数u,ud作为u和u'的解。
%在ode2solvertest版本中额外输入：真实解solution；输出：误差、迭代次数、区间数的信息。

%,diff,count,leafnum %,solution
%diff=zeros(1,30);
%count=0;
%leafnum=zeros(1,30);

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
u1=Func([a,c],0);
u2=Func([a,c],0);
u1d=Func([a,c],0);
u2d=Func([a,c],0);
nodeking=Node(a,c);
solveinv(nodeking,K,f,g1,g2,psi1,psi2);
evaluatenumber(nodeking,g1,g2);
nodeking.miul=0;
nodeking.miur=0;
solvesigma(nodeking);
Test=100;
leaflist=LeafList(nodeking);
while Test>Tol
    const=monitormax(nodeking)/(2^C);
    devide(nodeking,const,K,f,g1,g2,psi1,psi2);
    leaflist=LeafList(nodeking);
    upsweep(nodeking);
    downsweep(nodeking);
    l=max(size(leaflist));
    funcnode=zeros(1,l+1);
    funcnode(l+1)=c;
    coefu=zeros(K,l);
    coefud=zeros(K,l);
    J1=Jl(leaflist);
    J2=Jr(leaflist);
    for i=1:l
    solvesigma(leaflist(i));
    [w,wd]=localsolveequation(leaflist(i),K,J1(i),J2(i),ui,uid,s,g1,g2,g1d,g2d);
    coefu(:,i)=w';
    coefud(:,i)=wd';
    funcnode(i)=leaflist(i).a;
    end
    u1=u2;
    u1d=u2d;
    u2=Func(funcnode,coefu);
    u2d=Func(funcnode,coefud);
    Test=0;
    dif=@(t) (compute(u1,t)-compute(u2,t)).^2;
    for i=1:l
        Test=Test+integral(dif,funcnode(i),funcnode(i+1));
    end
    Test=sqrt(Test);
    summation=@(t) (compute(u1,t)+compute(u2,t)).^2;
    Test=Test/sqrt(integral(summation,a,c));

    %Test2=0;
    %dif2=@(t) (solution(t)-compute(u2,t)).^2;
    %for i=1:l
    %    Test2=Test2+integral(dif2,funcnode(i),funcnode(i+1));
    %end
    %count=count+1;
    %diff(count)=sqrt(Test2);
    %leafnum(count)=l;
end
u=u2;
ud=u2d;
end