function f=Jl(leaflist)
k=max(size(leaflist));
f=zeros(1,k);
f(1)=0;
for i=1:k-1
    f(i+1)=f(i)+leaflist(i).deltal+leaflist(i).miul*leaflist(i).alphal+leaflist(i).miur*leaflist(i).betal;
end
end