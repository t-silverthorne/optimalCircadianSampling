function r=get_perm_4tensor(N1,N2,N3,N4)
% Construction permutation tensor using recursion if N1*N2*N3*N4 is small
% enough, otherwise use nested for loops
if N1*N2*N3*N4<5e5
    r=randperm_recursive(NaN(N1*N3*N4,N2),1,N2,N1*N3*N4);
    r=reshape(r',N2,N1,N3,N4);
    r=pagetranspose(r);
else
    r=NaN(N1,N2,N3,N4);    
    for kk=1:N4
        for jj=1:N3
            for ii=1:N1
                r(ii,:,jj,kk)=randperm(N2,N2);
            end
        end
    end
end
end