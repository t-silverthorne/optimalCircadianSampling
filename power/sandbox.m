clf
clear
rng(3)
NL=4;
NR=4;
Nphi=2^4;
Nperm=1e3;
Amp=gpuArray(2.5);
Nbatches=5e2;
%nodeType='non-uniform';
simulatePWR(NL,NR,Nphi,Nperm,Amp,Nbatches,'uniform');
simulatePWR(NL,NR,Nphi,Nperm,Amp,Nbatches,'non-uniform');
function pwr=simulatePWR(NL,NR,Nphi,Nperm,Amp,Nbatches,nodeType)
Nmeas=NL+NR;
switch nodeType
    case 'uniform'
        t=linspace(0,1,Nmeas);
    case 'cheb'
        mc=1:Nmeas;
        t=cos((2*mc-1)*pi/2/Nmeas);
        t=(t+1)/2;
    case 'non-uniform'
        [~,t]=getSamplingSchedules(NL,NR,0,0.3);
end

for j=1:1
pwr=0;
for i=1:Nbatches

permMat=rand(Nphi,Nmeas,Nperm,'gpuArray');
[~,I]=sort(permMat,2);
% permMat=cell2mat(arrayfun(@(ii) cell2mat(arrayfun(@(ind) randperm(Nmeas),(1:Nphi)',"UniformOutput",false)),...
%     1:Nperm,...
%     'UniformOutput',false));
% I=gpuArray(reshape(permMat,[Nphi,Nmeas,Nperm]));

phi=gpuArray(transpose(linspace(0,2*pi,Nphi)));
eps=randn(Nphi,Nmeas,'gpuArray');

A1=Amp*sin(phi);
A2=Amp*cos(phi);

Y=A1.*sin(2*pi*t)+A2.*cos(2*pi*t)+eps;
SSres_obs=getSSres(Y,t);

zts=t;
x1=gpuArray(sin(2*pi*t));
x2=gpuArray(cos(2*pi*t));
x0=gpuArray(ones(1,size(Y,2)));
X= [x0' x1' x2'];

m=size(Y,1);n=size(Y,2);
offMat=repmat((0:m-1)',1,n)*n;
Yp=Y';
YI=pagetranspose(Yp(pagetranspose(I+offMat)));
betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr=pwr+sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nphi;
end
pwr=pwr/Nbatches;
fprintf('%f\n',pwr)
end
% count=0;
% for ii=1:Nperm
%     [~,I]=sort(permMat(:,:,ii),2);
%     count=count+(getSSres(Y(I),t)>SSres_obs);
% end
% 
% transpose(count/Nperm)
end


function SSres=getSSres(Y,t)
zts=t;
x1=gpuArray(sin(2*pi*t));
x2=gpuArray(cos(2*pi*t));
x0=gpuArray(ones(1,size(Y,2)));
X= [x0' x1' x2'];

betas=(X'*X)\(X'*Y');

fits=(X*betas)';
% size(fits)
% plot(fits(1,:))
% hold on
% plot(Y(1,:))

SSres=sqrt(sum((fits-Y).^2,2));
end
