function [u,ud,diff,count]=nonlinearode2solver2test(a,c,f,fu,fud,nodelist,K,B,Cons,Tol,epsilon,solution)
count=0;
diff=zeros(1,100);

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
Test2=0;
dif2=@(t) (solution(t)-compute(u,t)).^2;
for i=1:l
    Test2=Test2+integral(dif2,nodelist(i),nodelist(i+1));
end
diff(count)=sqrt(Test2);

end
end