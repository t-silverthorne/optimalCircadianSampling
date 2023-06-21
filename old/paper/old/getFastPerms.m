function YI = getFastPerms(p,Y)
d=p.Nmeas;
R=NaN(p.Nperm,d);
for ii=1:p.Nperm
    R(ii,:)=randperm(d,d);
end
R;
N3=p.Nperm;
N4=p.Nacro;
%max(sum(pagetranspose(reshape(R',[d 1 N3]))))
%min(sum(pagetranspose(reshape(R',[d 1 N3]))))

R=pagetranspose(reshape(R',[d 1 p.Nperm]));
I=eye(d);

ind3=repmat(ones(d,d),[1 1 N3 N4]);
for ii=1:N3
    ind3(:,:,ii,:)=ind3(:,:,ii,:)*ii;
end
ind4=repmat(ones(d,d),[1 1 N3 N4]);
for ii=1:N4
    ind4(:,:,:,ii)=ind4(:,:,:,ii)*ii;
end
I=repmat(I,[1,1,p.Nperm,p.Nacro]);
YI=pagemtimes(Y, I(sub2ind(size(I),repmat((1:d)',[1 d N3 N4]),repmat(R,[d 1 1 N4]),ind3,ind4)) );

end

