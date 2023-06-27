format long

epsilon=0.001;
v=100;
p=@(t) 0;
q=@(t) 0;
f=@(t) 2;
K=10;
a=-1;
c=1;
solution=@(t) 0;
A=[1,1;0,0];
gamma1=-1;
gamma2=1;
C=4;
Tol=10^(-9);
[u,ud,diff,count,leafnum]=ode2solvertest(a,c,K,p,q,f,A,gamma1,gamma2,C,Tol,solution);

%N=8;
%u;
%ud;
%diff=zeros(1,N);
%for i=1:N
%nodelist=a:(c-a)/(2^i):c;
%[u,ud,diff(i)]=directode2solvertest(a,c,K,p,q,f,A,gamma1,gamma2,nodelist,solution);
%end

xx=a:0.01:c;
num=max(size(xx));
yy=zeros(1,num);
for i=1:num
    yy(i)=compute(u,xx(i));
end
figure(1)
plot(xx,solution(xx));
figure(2)
plot(xx,solution(xx),'red',xx,yy,'green');
max(abs(solution(xx)-yy))

count
diff
figure(3)
semilogy(1:count,diff(1:count));
figure(4)
plot(1:count,leafnum(1:count));

%figure(3)
%semilogy(1:N,diff);