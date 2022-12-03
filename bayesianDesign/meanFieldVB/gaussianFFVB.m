% assume posterior is sum of mutlivariate normal distributions
% syms d N
% simplify(solve(N-d*(d+1)/2,d))
% f=@(N) (8*N + 1)^(1/2)/2 - 1/2
% f()
clear 
% GPU random variable generation
NS=1e7;% number of samples
d=3;
Ngauss=3;

Lvec = rand(d*(d+1)/2,Ngauss);
muvec = rand(d,Ngauss);
wvec = rand(Ngauss,1); wvec = wvec/sum(wvec);
fprintf('running\n')
gputimeit(@() generateSamples(muvec,Lvec,wvec,NS,d,Ngauss,true))
timeit(@() generateSamples(muvec,Lvec,wvec,NS,d,Ngauss,false))
function S=generateSamples(muvec,Lvec,wvec,NS,d,Ngauss,useGPU)

switch useGPU
    case true
        eps  = randn(d,1,NS,'gpuArray');
        u    = rand(1,1,NS,'gpuArray');
    case false
        eps  = randn(d,1,NS);
        u    = rand(1,1,NS);
end

avec = cumsum(wvec);
avecm1= [0; avec(1:end-1)];
xi=avecm1<u & u<avec;
counts=sum(xi,2);
S=NaN(d,NS,'gpuArray');
cc=cumsum(counts);
Lmatnow=ivh(Lvec(:,1));
S(:,1:cc(1))=muvec(:,1) + Lmatnow*eps(:,1:cc(1));
for ii=2:Ngauss
    i0=cc(ii-1)+1;
    Lmatnow=ivh(Lvec(:,ii));
    S(:,i0:cc(ii))=muvec(:,ii) + Lmatnow*eps(:,i0:cc(ii));
end
end
   
function Avech=vech(A)
d=size(A,1);
Avech=[];
for ii=1:size(A,2)
    for jj=ii:d
        Avech(end+1)=A(jj,ii);
    end
end
end


function Lmat = ivh(Lvec)
d= (8*length(Lvec) + 1)^(1/2)/2 - 1/2;
Lmat = zeros(d,d);
i0=1;
for i=1:d
    Lmat(i:d,i)=Lvec(i0:i0+d-i);
    i0=i0+d+1-i;
end
end

function S=sample_from_q(Nsamp,muCell,cholMatCell,wvec)
% input:
% muCell    - list of means
% sigmaCell - list of covariance matrices
% wvec      - weights for each gaussian 
% output: samples from the summed distribution
% generate random samples from qlambda(theta)
mu_inds = cumsum(mnrnd(Nsamp,wvec));
S=randn([Nsamp,3]);

for ii=1:length(mu_inds)
    R=cholMatCell{ii};
    if ii>1
        Sstart=mu_inds(ii-1)+1;
    else
        Sstart=1;
    end
    Sstop=mu_inds(ii);
    S(Sstart:Sstop)=muCell{ii} + S(Sstart:Sstop)*R;
end
end