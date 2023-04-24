p.Nmeas=10;
p.Nacro=32;
p.Nresidual=1e3;
p.Nperm=1e2;
p.freq=1;
p.Amp=1;
p.noise=.5;

[t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
R=getPermutations(p.Nresidual,p.Nmeas,p.Nperm,p.Nacro);
Y=getSimulatedData(t,p);

[I3,I4]=constructUtilMats(p);

tic
R=getPermutations(p.Nresidual,p.Nmeas,p.Nperm,p.Nacro);
getPermutedData(Y,R,I3,I4);
toc

tic
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
toc
%%
size(repmat((1:d)',[1 d N3 N4]))
size(repmat(R,[d 1 1 N4]))
size(ind3)
size(ind4)