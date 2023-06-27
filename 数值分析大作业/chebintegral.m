function f=chebintegral(K,sigma,b1,b2)
f=0;
for i=1:2:K
    f=f+sigma(i)*2/(1-(i-1)^2);
end
f=f*(b2-b1)/2;
end