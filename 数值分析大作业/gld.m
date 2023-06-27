function f=gld(x,A,a,c)
q=Q0(A);
if q==-1
    f=A(2,1)*sinh(x-a)-A(1,1)*cosh(x-a);
else
    f=A(1,1);
end
end