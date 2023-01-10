% test if uniform vs non-uniform sampling improves performance when data is
% generated using randomly obtained Fourier coefficients

clf
clear
param.NL=5;
param.NR=3;
param.Nperm=1e3;
param.Nbatches=5e2; % SMALL RIGHT NOW
param.per=1; % period used in regression model

param.NfourierSamps=100; % num. fourier samples
param.fourierMeans=[5 3.1]*4;
param.fourierSigma=[1 1];
param.NfourierComps=length(param.fourierMeans);

simulatePWR(param,'uniform');
simulatePWR(param,'non-uniform');

function pwr=simulatePWR(param,nodeType)
NL=param.NL;
NR=param.NR;
Nperm=param.Nperm;
Nbatches=param.Nbatches;

NfourierComps=param.NfourierComps;
NfourierSamps=param.NfourierSamps;
fourierMeans=param.fourierMeans;
fourierSigma=param.fourierSigma;

Nmeas=NL+NR;
per=param.per;
switch nodeType
    case 'uniform'
        t=linspace(0,1,Nmeas);
    case 'cheb'
        mc=1:Nmeas;
        t=cos((2*mc-1)*pi/2/Nmeas);
        t=(t+1)/2;
    case 'non-uniform'
        r=0;%rand*0.7;
        [~,t]=getSamplingSchedules(NL,NR,r,r+1/3);
        t=gpuArray(t);
end

pwr=0;
for i=1:Nbatches

permMat=rand(NfourierSamps,Nmeas,Nperm,'gpuArray');
[~,I]=sort(permMat,2);

eps=randn(NfourierSamps,Nmeas,'gpuArray');

fourierSigma=reshape(fourierSigma,1,1,length(fourierSigma));
fourierMeans=reshape(fourierMeans,1,1,length(fourierMeans));
theta=2*pi*rand(NfourierSamps,1,NfourierComps,'gpuArray');
ampvec=fourierMeans+fourierSigma.*randn(NfourierSamps,1,NfourierComps,'gpuArray');

pervec=1:NfourierComps;
pervec=reshape(pervec,1,1,length(pervec));
Y=sum(ampvec.*cos(2*pi*t.*pervec-theta),3)+eps;

SSres_obs=getSSres(Y,per,t);

x1=gpuArray(sin(2*pi*per*t));
x2=gpuArray(cos(2*pi*per*t));
x0=gpuArray(ones(1,size(Y,2)));
X= [x0' x1' x2'];

m=size(Y,1);n=size(Y,2);
offMat=repmat((0:m-1)',1,n)*n;
Yp=Y';
YI=pagetranspose(Yp(pagetranspose(I+offMat)));
betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr=pwr+sum(sum(SSres>SSres_obs,3)/Nperm>.95)/NfourierSamps;
end
pwr=pwr/Nbatches;
fprintf('%f\n',pwr)
end



function SSres=getSSres(Y,per,t)

x1=gpuArray(sin(2*pi*per*t));
x2=gpuArray(cos(2*pi*per*t));
x0=gpuArray(ones(1,size(Y,2)));
X= [x0' x1' x2'];

betas=(X'*X)\(X'*Y');

fits=(X*betas)';

SSres=sqrt(sum((fits-Y).^2,2));
end

% param.NL=20;
% param.NR=10;
% param.Nperm=1e3;
% param.Amp=gpuArray(2.5);
% param.Nbatches=1e2;
% param.per=1; % period used in regression model
% 
% param.NfourierComps=3; % num. fourier components included in data generation
% param.NfourierSamps=100; % num. fourier samples
% param.fourierMeans=gpuArray([1 1 1]*2);
% param.fourierSigma=gpuArray([1 1 1]);
% 
% simulatePWR(param,'uniform');
% simulatePWR(param,'non-uniform');
