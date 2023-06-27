    function del=Delta(t,omiga,phidelta1,phidelta1d,g1,g)
        del=zeros(2,max(size(t)));
        for i=1:max(size(t))
        del(:,i)=-omiga(t(i))*[compute(phidelta1,t(i));compute(phidelta1d,t(i))-g1(t(i))]+g(t(i));
        end
    end