function f=cheb(a,c,K)
f=(c-a)*(cos(((2*K-1):-2:1)*pi/(2*K)))/2+(c+a)/2;
end