p=@(t) 2./t;
q=@(t) 1./(t.^4);
f=@(t) 0;
solution=@(t) cos(1./t);
a=0.02;
c=1;
A=[1,1;0,0];
gamma1=cos(1/a);
gamma2=cos(1/c);
K=10;
C=4;
Tol=10^(-40);
[u,ud,diff,count,leafnum]=ode2solvertest(a,c,K,p,q,f,A,gamma1,gamma2,C,Tol,solution);

%N=8;
%u;
%ud;
%diff=zeros(1,N);
%for i=1:N
%nodelist=a:(c-a)/(2^i):c;
%[u,ud,diff(i)]=directode2solvertest(a,c,K,p,q,f,A,gamma1,gamma2,nodelist,solution);
%end

xx=a:0.001:c;
num=max(size(xx));
yy=zeros(1,num);
for i=1:num
    yy(i)=compute(u,xx(i));
end
figure(1)
plot(xx,yy)
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

