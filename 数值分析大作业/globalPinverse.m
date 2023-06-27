function sigma=globalPinverse(f,K,nodelist,g1,g2,psi1,psi2)
l=max(size(nodelist))-1;
A=zeros(K*l);
value=zeros(K*l,1);
for node=1:l
    t=cheb(nodelist(node),nodelist(node+1),K);
    for i=1:K
        value((node-1)*K+i)=f(t(i));
    end
end
for nodej=1:l
    b1=nodelist(nodej);
    b2=nodelist(nodej+1);
    for j=1:K
        prodl=@(t) g1(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
        prodr=@(t) g2(t).*cos((j-1)*acos(2*(t-b1)/(b2-b1)-1));
        prodlint=integral(prodl,b1,b2);
        prodrint=integral(prodr,b1,b2);
        c1=intl(prodl,K,b1,b2);
        c2=intl(prodr,K,b1,b2);
        for nodei=1:l
            t=cheb(nodelist(nodei),nodelist(nodei+1),K);
            for i=1:K
                if nodei<nodej
                    A((nodei-1)*K+i,(nodej-1)*K+j)=psi2(t(i))*prodrint;
                else 
                if nodei==nodej
                    A((nodei-1)*K+i,(nodej-1)*K+j)=cos((j-1)*acos(2*(t(i)-b1)/(b2-b1)-1))+psi1(t(i))*chebvalue(c1,b1,b2,t(i))+psi2(t(i))*chebvalue(c2,b1,b2,t(i));
                else 
                    A((nodei-1)*K+i,(nodej-1)*K+j)=psi1(t(i))*prodlint;
                end
                end

            end
        end
    end
end
w=A\value;
sigma=zeros(K,l);
for i=1:l
    sigma(:,i)=w(K*(i-1)+1:K*i,1);
end
end