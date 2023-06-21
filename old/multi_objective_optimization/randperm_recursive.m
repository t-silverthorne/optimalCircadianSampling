function R=randperm_recursive(R,ind,d,N)
R(ind,:)=randperm(d,d);
if ind<N
    R=randperm_recursive(R,ind+1,d,N);
end 
end