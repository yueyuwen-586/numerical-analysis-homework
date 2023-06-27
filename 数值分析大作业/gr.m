function f=gr(x,A,a,c)
q=Q0(A);
if q==-1
    f=A(2,2)*cosh(x-c)-A(1,2)*sinh(x-c);
else
    f=A(1,2)*(x-c)-A(2,2);
end
end