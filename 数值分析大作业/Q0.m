function f=Q0(A)
if abs(A(1,1))<abs(A(2,1)) && abs(A(1,2))<abs(A(2,2))
    f=-1;
else
    f=0;
end
end