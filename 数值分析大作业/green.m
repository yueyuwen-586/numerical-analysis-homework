function w=green(A,C,fund,a,c,x,t)
%计算格林函数
D=A*fund(a)+C*fund(c);
s=size(fund(a),1);
I=eye(s);
J=@(u) -(I/D)*C*fund(c)*(I/fund(u));
if t<x
w=fund(x)*(I/fund(t)+J(t));
else
w=fund(x)*J(t);
end