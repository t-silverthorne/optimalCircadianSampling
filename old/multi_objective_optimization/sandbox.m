
% 2 pages, use same permutation on both
Y=rand(3,3,2);

r=[randperm(3,3);randperm(3,3);randperm(3,3)];

Y

ii(:,:,1)=ones(3,3);
ii(:,:,2)=2*ones(3,3);
Y(sub2ind(size(Y),repmat([1;2;3],1,3,2),repmat(r,1,1,2),ii))
%%
N1=1e3;
N2=8;
N3=1e2;
N4=30;
Y=rand(N1,N2,N3,N4);

tic
r=get_perm_tensor(N1,N2,N3,N4);

I3=NaN(N1,N2,N3);
I4=NaN(N1,N2,N3,N4);
for ii=1:N3
    I3(:,:,ii)=ii*ones(N1,N2);
end
for ii=1:N4
    I4(:,:,:,ii)=ii*ones(N1,N2,N3);
end

YI=Y(sub2ind(size(Y),repmat((1:N1)',1,N2,N3,N4),r,repmat(I3,1,1,1,N4),I4));
toc

function YI=get_permuted_Y_rank4(Y)
N1=size(Y,1);N2=size(Y,2);N3=size(Y,3);N4=size(Y,4);
r=get_perm_tensor(N1,N2,N3,N4);
I3=NaN(N1,N2,N3);
I4=NaN(N1,N2,N3,N4);
for ii=1:N3
    I3(:,:,ii)=ii*ones(N1,N2);
end
for ii=1:N4
    I4(:,:,:,ii)=ii*ones(N1,N2,N3);
end
YI=Y(sub2ind(size(Y),repmat((1:N1)',1,N2,N3,N4),r,repmat(I3,1,1,1,N4),I4));
end

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