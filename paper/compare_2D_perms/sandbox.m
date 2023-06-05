N=5e5;
d=8;
% vectorized FY method
Ydat=rand(N,d);
useGPU=false;

if useGPU
    Ydat=gpuArray(Ydat);
end
Ydat_init=Ydat;
tic
for ii=1:d-1
    inds=randsample((ii:d)',N,true);
    switch useGPU
        case true
            A=gpuArray(zeros(N,d));
            B=gpuArray(ones(N,d));  % everything not involved in permutation
            I=gpuArray(eye(d,d));
        case false
            A=zeros(N,d);
            B=ones(N,d);  % everything not involved in permutation
            I=eye(d,d);
    end
    A(sub2ind(size(A),(1:N)',inds))=1;
    B(:,ii)=0;
    B(sub2ind(size(B),(1:N)',inds))=0;
    A(:,ii)=A(:,ii)/2; % avoid double counting
    
    %Ydat
    Ydat=Ydat.*B +  A.*Ydat(:,ii)  + ( ones(N,1)*I(ii,:) ).*( sum(Ydat.*A,2) );
end
toc
Ydat=gather(Ydat);
tic
for ii=1:N
    Ydat(ii,:) =Ydat(ii,randperm(d,d));
end
toc