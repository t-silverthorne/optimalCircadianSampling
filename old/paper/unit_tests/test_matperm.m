
addpath('utils_core')
N3=2;d=5;N4=2;
R=getPermutations(1,d,N3,N4);
R;
I=eye(d);
I=repmat(I,[1 1 N3 N4]);

ind3=repmat(ones(d,d),[1 1 N3 N4]);
for ii=1:N3
    ind3(:,:,ii,:)=ind3(:,:,ii,:)*ii;
end
ind4=repmat(ones(d,d),[1 1 N3 N4]);
for ii=1:N4
    ind4(:,:,:,ii)=ind4(:,:,:,ii)*ii;
end

I(sub2ind(size(I),repmat((1:d)',[1 d N3 N4]),repmat(R,[d 1 1 1]),ind3,ind4))

%%
size(repmat((1:d)',[1 d N3]))
size(R)
size(ind3)
%I(sub2ind(R(1:2,:))