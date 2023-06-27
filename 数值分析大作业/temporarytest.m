%temporary test
help trans
%trans([0,0.5,1],@(t) t,10)

funcnode=[0,0.5,1];
f=@(t) t;
K=10;
l=max(size(funcnode))-1;
            b1=funcnode(1);
            b2=funcnode(l);
            A=zeros(K,l);
            for i=1:l
                v=f(cheb(funcnode(i),funcnode(i+1),K));
                A(:,i)=(chebint(v,b1,b2,K))';
            end
            func=Func(funcnode,A)