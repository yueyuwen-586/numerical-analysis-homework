function [u,ud]=nonlinearode2solver(a,c,f,fu,fud,m,K,B,Cons,Tol,epsilon)
A=[B(1,1),B(2,1);0,0];
C=[0,0;B(1,2),B(2,2)];
Dn=A+C;
if rank(Dn)==2
p0=@(t) 0;
fund=@(t) eye(2);
G0=@(x,t) green(A,C,fund,a,c,x,t);
elseif sin(a-c)~=0
p0=@(t) [0,-1;1,0];
fund=@(t) [cos(t),-sin(t);sin(t),cos(t)];
G0=@(x,t) green(A,C,fund,a,c,x,t);
else 
p0=@(t) [0,-1;2,0];
fund=@(t) [cos(sqrt(2)*t),-sin(sqrt(2)*t);sin(sqrt(2)*t),cos(sqrt(2)*t)];
G0=@(x,t) green(A,C,fund,a,c,x,t);
end
funcnode=a:(c-a)/(2^m):c;
deltaL2=100;
sigmaL2=1;

Phi1=Func([a,c],0);
Phi2=Func([a,c],0);
Sigma1=Func([a,c],0);
Sigma2=Func([a,c],0);
while deltaL2/sigmaL2>epsilon
    %phi=Phi
    phi1=@(t) compute(Phi1,t);
    phi2=@(t) compute(Phi2,t);
    %sigma=Sigma
    sigma1=@(t) compute(Sigma1,t);
    sigma2=@(t) compute(Sigma2,t);
    sigma=@(t) [sigma1(t);sigma2(t)];
    %omiga
    fphi=@(t) f(t,phi1(t),phi2(t));
    fuphi=@(t) fu(t,phi1(t),phi2(t));
    fudphi=@(t) fud(t,phi1(t),phi2(t));
    omiga=@(t) [0,-1;fuphi(t),fudphi(t)]-p0(t);
    %g
    g=@(t) p0(t)*[phi1(t);phi2(t)]+[phi2(t);fphi(t)]-sigma(t);
    g1=@(t) takeone(g,1,t);
    g2=@(t) takeone(g,2,t);
    h=10^(-5);
    g1d=@(t) (g1(t+h)-g1(t-h))/(2*h);
    %delta(use phidelta)
    p=@(t) -fudphi(t);
    q=@(t) -fuphi(t);
    f0=@(t) g2(t)+g1d(t)-fuphi(t).*g1(t);
    Boundry=B;
    gamma1=B(2,1)*([1,0]*g(a));
    gamma2=B(2,2)*([1,0]*g(c));
    [phidelta1,phidelta1d]=ode2solver(a,c,K,p,q,f0,Boundry,gamma1,gamma2,Cons,Tol);
    delta=@(t) Delta(t,omiga,phidelta1,phidelta1d,g1,g);
    %sigma=sigma+delta
    sum1=@(t) compute(Sigma1,t)+takeone(delta,1,t);
    sum2=@(t) compute(Sigma2,t)+takeone(delta,2,t);
    Sigma1=Func(funcnode,generatecoef(funcnode,sum1,K));
    Sigma2=Func(funcnode,generatecoef(funcnode,sum2,K));
    %Sigma1=Sigma1copy;
    %Sigma2=Sigma2copy;

    %phi
    prod=@(x,t) G0(x,t)*[compute(Sigma1,t);compute(Sigma2,t)];
    phi1v=@(x) integral(@(t) takeone(prod(x,t),1),a,c);
    phi2v=@(x) integral(@(t) takeone(prod(x,t),2),a,c);
    Phi1=Func(funcnode,generatecoef(funcnode,phi1v,K));
    Phi2=Func(funcnode,generatecoef(funcnode,phi2v,K));
    %computeL2
    sqr1=@(t) compute(Sigma1,t)^2+compute(Sigma2,t)^2;
    sqr2=@(t) delta(t)'*delta(t);
    deltaL2=integral(sqr1,a,c);
    deltaL2=sqrt(deltaL2);
    sigmaL2=integral(sqr2,a,c);
    sigmaL2=sqrt(sigmaL2);
end

u=Phi1;
ud=Phi2;
end