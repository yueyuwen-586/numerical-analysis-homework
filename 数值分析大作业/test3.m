f=@(t,u,ud) 2*(u-t/2+1.5).^3;
fu=@(t,u,ud) 6*(u-t/2+1.5).^2;
fud=@(t,u,ud) 0;
K=10;
a=1;
c=2;
B=[1,1;0,0];
solution=@(t) 1./t+t/2-3/2;
Cons=4;
Tol=10^(-10);
m=5;
nodelist=a:(c-a)/(2^m):c;
epsilon=10^(-3);
[u,ud,diff,count]=nonlinearode2solver2test(a,c,f,fu,fud,nodelist,K,B,Cons,Tol,epsilon,solution);

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
figure(3)
semilogy(1:count,diff(1:count));