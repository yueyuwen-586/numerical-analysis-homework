function w=chebvalue(c,b1,b2,x)
k=max(size(c));
w=c*cos((0:k-1)'*acos(2*(x-b1)/(b2-b1)-1));
end