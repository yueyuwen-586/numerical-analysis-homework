function f=psil(x,p,q,A,a,c)
f=(p(x).*grd(x,A,a,c)+q(x).*gr(x,A,a,c))/S(A,a,c);
end