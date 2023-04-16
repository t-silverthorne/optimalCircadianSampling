function Ypages = get_permuted_pages_recursive(Ypages,ind,fast)
N=size(Ypages,1);
d=size(Ypages,2);
Nperm=size(Ypages,3);

if fast
    R=randperm_recursive(NaN(1,d),1,d,1);
    R=repmat(R,N,1);
else
    R=randperm_recursive(NaN(N,d),1,d,N);
end
Ypages(:,:,ind)=Ypages(sub2ind(size(R),repmat((1:N)',1,d),R));
if ind<Nperm
    Ypages=get_permuted_pages_recursive(Ypages,ind+1,fast);
end
end