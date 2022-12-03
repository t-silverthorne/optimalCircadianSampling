% assume posterior is sum of mutlivariate normal distributions
% syms d N
% simplify(solve(N-d*(d+1)/2,d))
% f=@(N) (8*N + 1)^(1/2)/2 - 1/2
% f()
clear all
% GPU random variable generation
NS=1e7;% number of samples
d=3;
Ngauss=3;

Lvec = rand(d*(d+1)/2,Ngauss);
muvec = rand(d,Ngauss);
wvec = rand(Ngauss,1); wvec = wvec/sum(wvec);

tic
isgpuarray(generateSamples(muvec,Lvec,wvec,NS,d,Ngauss))
toc

function S=generateSamples(muvec,Lvec,wvec,NS,d,Ngauss)
useGPU=true;
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
xi=avecm1<u & u<avec ;

S=(muvec(:,1)+pagemtimes(ivh(Lvec(:,1)),eps) ).*xi(1,:,:);

if Ngauss>1
    S=S+(muvec(:,2)+pagemtimes(ivh(Lvec(:,2)),eps) ).*xi(2,:,:);
end
if Ngauss>2
    S=S+(muvec(:,3)+pagemtimes(ivh(Lvec(:,3)),eps) ).*xi(3,:,:);
end
if Ngauss>3
    S=S+(muvec(:,4)+pagemtimes(ivh(Lvec(:,4)),eps) ).*xi(4,:,:);
end
if Ngauss>4
    S=S+(muvec(:,5)+pagemtimes(ivh(Lvec(:,5)),eps) ).*xi(5,:,:);
end
if Ngauss>5
    S=S+(muvec(:,6)+pagemtimes(ivh(Lvec(:,6)),eps) ).*xi(6,:,:);
end
if Ngauss >6
    print("stop: too many gaussians")
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