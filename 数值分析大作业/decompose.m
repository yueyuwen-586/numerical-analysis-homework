function f=decompose(A,a,c,gamma1,gamma2)
B=[a*A(1,1)+A(2,1),A(1,1);c*A(1,2)+A(2,2),A(1,2)];
g=B\[gamma1;gamma2];
f=g';
end