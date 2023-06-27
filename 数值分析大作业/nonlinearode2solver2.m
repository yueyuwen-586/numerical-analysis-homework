function [u,ud]=nonlinearode2solver2(a,c,f,fu,fud,nodelist,K,B,Cons,Tol,epsilon)
%这是非线性方程求解器，输入区间端点a、c，插值次数K，插值区间端点,误差允许范围epsilon
%二阶方程参数函数f，偏导f_u,f_u'，边值条件矩阵B=[zetal0,zetar0;zetal1,zetar1]，
%线性算法中的C和Tol。输出Func类函数u,ud作为u和u'的解。
%在nonlinearode2solver2test版本中额外输入：真实解solution；输出：误差、迭代次数的信息。

count=0;
l=max(size(nodelist))-1;
u=Func(nodelist,zeros(K,l));
ud=Func(nodelist,zeros(K,l));
udd=Func(nodelist,zeros(K,l));
vL2=100;
uL2=1;
while vL2/uL2>epsilon
p=@(t) -fud(t,compute(u,t),compute(ud,t));
q=@(t) -fu(t,compute(u,t),compute(ud,t));
f0=@(t) -compute(udd,t)+f(t,compute(u,t),compute(ud,t));
[v,vd]=ode2solver(a,c,K,p,q,f0,B,0,0,Cons,Tol);
vddfunc=@(t) -compute(vd,t).*p(t)-compute(v,t).*q(t)+f0(t);
ufunc=@(t) compute(u,t)+compute(v,t);
udfunc=@(t) compute(ud,t)+compute(vd,t);
uddfunc=@(t) compute(udd,t)+vddfunc(t);
coefu=generatecoef(nodelist,ufunc,K);
coefud=generatecoef(nodelist,udfunc,K);
coefudd=generatecoef(nodelist,uddfunc,K);
u=Func(nodelist,coefu);
ud=Func(nodelist,coefud);
udd=Func(nodelist,coefudd);
sqr1=@(t) compute(v,t).^2;
sqr2=@(t) compute(u,t).^2;
vL2=sqrt(integral(sqr1,a,c));
uL2=sqrt(integral(sqr2,a,c));
count=count+1;
end
end