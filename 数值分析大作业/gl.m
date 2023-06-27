function f=gl(x,A,a,c)
q=Q0(A);
if q==-1
    f=A(2,1)*cosh(x-a)-A(1,1)*sinh(x-a);
else
    f=A(1,1)*(x-a)-A(2,1);
end
end