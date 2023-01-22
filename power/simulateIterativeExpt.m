% test if uniform vs non-uniform sampling improves performance when data is
% generated using randomly obtained Fourier coefficients
close all
clear
param.NL=5;
param.NR=3;
simtype='fast';
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'long'
        param.Nperm=1e3;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'verylong'
        param.Nperm=1e3;
        param.Nresidual=5e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
end

% Uniform oscillates, others do not
% param.NL=8; param.NR=4;
% param.per=18; % period used in regression model
% param.pertrue=param.per;
% param.Amp=1.5;
Numdelta=1;
tiledlayout(1,Numdelta)
deltavals=0;%linspace(-5,5,Numdelta);
for j=1:Numdelta
    nexttile(j)
    param.per=2; % period used in regression model
    param.pertrue=(1+deltavals(j)/100)*param.per;
    param.Amp=2;
    
    pwrUnif=simulatePWR(param,'uniform');
    plot(linspace(0,2*pi,param.Nacro),pwrUnif)
    hold on
    
    pwrNU=simulatePWR(param,'non-uniform');
    plot(linspace(0,2*pi,param.Nacro),pwrNU)
    
    pwrCH=simulatePWR(param,'cheb');
    plot(linspace(0,2*pi,param.Nacro),pwrCH)
    legend({'unif','nu','cheb'})
    ylim([0,1])
    drawnow
end
function pwr=simulatePWR(param,nodeType)
NL=param.NL;
NR=param.NR;
Nperm=param.Nperm;
Nresidual=param.Nresidual;

Nacro=param.Nacro;
Amp=param.Amp;
Nmeas=NL+NR;
per=param.per;
pertrue=param.pertrue;

switch nodeType
    case 'uniform'
        [t,~]=getSamplingSchedules(NL,NR,0,0.5);
        %xline(t*2*pi,'Color',[1 0 0])
        %hold on
    case 'cheb'
        mc=1:Nmeas;
        t=cos((2*mc-1)*pi/2/Nmeas);
        t=(t+1)/2;
        t=sort(t);
        %xline(t*2*pi,'Color',[0 1 0])
        %hold on
    case 'non-uniform'
        [~,t]=getSamplingSchedules(NL,NR,0,0.5);
        %xline(t*2*pi,'Color',[0 0 1])
        %hold on
end

hold on
pwr=[];
acrovec=linspace(0,2*pi,Nacro);
for i=1:Nacro
acro=acrovec(i);

permMat=rand(Nresidual,Nmeas,Nperm,'gpuArray');
[~,I]=sort(permMat,2);

eps=randn(Nresidual,Nmeas,'gpuArray');
Y=Amp*cos(2*pi*t*pertrue-acro)+eps;

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
pwr(end+1)=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nresidual;
end
%pwr=pwr/Nresidual;
%fprintf('%f\n',pwr)
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
% param.Nresidual=1e2;
% param.per=1; % period used in regression model
% 
% param.NfourierComps=3; % num. fourier components included in data generation
% param.Nacro=100; % num. fourier samples
% param.fourierMeans=gpuArray([1 1 1]*2);
% param.fourierSigma=gpuArray([1 1 1]);
% 
% simulatePWR(param,'uniform');
% simulatePWR(param,'non-uniform');
