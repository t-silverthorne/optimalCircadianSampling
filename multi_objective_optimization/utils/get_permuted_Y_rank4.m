function YI=get_permuted_Y_rank4(Y)
N1=size(Y,1);N2=size(Y,2);N3=size(Y,3);N4=size(Y,4);
r=get_perm_4tensor(N1,N2,N3,N4);
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
