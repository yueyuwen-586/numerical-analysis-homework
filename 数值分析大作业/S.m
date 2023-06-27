function f=S(A,a,c)
if Q0(A)==0
    f=det(A)+A(1,1)*A(1,2)*(c-a);
else
    f=det(A)*cosh(a-c)+(A(2,1)*A(2,2)-A(1,1)*A(1,2))*sinh(a-c);
end
end