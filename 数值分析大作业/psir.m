function f=psir(x,p,q,A,a,c)
f=(p(x).*gld(x,A,a,c)+q(x).*gl(x,A,a,c))/S(A,a,c);
end