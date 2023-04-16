function Pmat = construct_all_perms(d)
% Pmat is a rank 3 array of all permutation matrices on d elements,
% i.e. each page of Pmat is a permutation matrix
pall=perms(1:d);
Pmat=NaN(d,d,factorial(d));
I=eye(d);
for ii=1:factorial(d)
    Pmat(:,:,ii)= I(:,pall(ii,:)); % might be possible to vectorize
end
end

