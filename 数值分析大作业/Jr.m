function f=Jr(leaflist)
k=max(size(leaflist));
f=zeros(1,k);
f(k)=0;
for i=k:-1:2
    f(i-1)=f(i)+leaflist(i).deltar+leaflist(i).miul*leaflist(i).alphar+leaflist(i).miur*leaflist(i).betar;
end
end