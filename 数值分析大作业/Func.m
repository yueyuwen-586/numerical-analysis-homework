classdef Func<handle
    properties
        a %节点序列
        Coef %切比雪夫多项式系数

        %用分段插值来记录函数，在[a(i),a(i+1)]上，切比雪夫插值多项式
        %系数向量为矩阵Coef的第i列。由compute函数计算其在x处的值。
    end

    methods
        function func=trans(funcnode,f,K)
            l=max(size(funcnode))-1;
            b1=funcnode(1);
            b2=funcnode(l);
            A=zeros(K,l);
            for i=1:l
                v=f(cheb(funcnode(i),funcnode(i+1),K));
                A(:,i)=(chebint(v,b1,b2,K))';
            end
            func=Func(funcnode,A);
        end
      function f=compute(func,x)
          i=1;
          while (func.a(i+1))<x 
              i=i+1;
          end
          b1=func.a(i);
          b2=func.a(i+1);
          coefi=func.Coef(:,i)';
          f=chebvalue(coefi,b1,b2,x);
      end
      function func = Func(a,coef)
         if (nargin > 0)
            func.a=a;
            func.Coef=coef;
         end
      end
    end
    methods (Access = private)
      function delete(~)
      end
    end
end