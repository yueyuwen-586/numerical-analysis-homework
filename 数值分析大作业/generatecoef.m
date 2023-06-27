function A=generatecoef(funcnode,f,K)
            l=max(size(funcnode))-1;
            b1=funcnode(1);
            b2=funcnode(l);
            A=zeros(K,l);
            for i=1:l
                v=f(cheb(funcnode(i),funcnode(i+1),K));
                A(:,i)=(chebint(v,b1,b2,K))';
            end
end