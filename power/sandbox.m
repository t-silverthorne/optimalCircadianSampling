clf
clear
rng(3)
Nmeas=8;
Nphi=2^4;
Nperm=1e4;
Amp=gpuArray(2.6);
Nbatches=1e2;
nodeType='uniform'
for j=1:5
pwr=0;
for i=1:Nbatches
permMat=rand(Nphi,Nmeas,Nperm,'gpuArray');

%phi=2*pi*rand(Nphi,1,'gpuArray');
phi=gpuArray(transpose(linspace(0,2*pi,Nphi)));
eps=randn(Nphi,Nmeas,'gpuArray');

A1=Amp*sin(phi);
A2=Amp*cos(phi);
switch nodeType
    case 'uniform'
        t=linspace(0,1,Nmeas);
    case 'cheb'
        m=1:Nmeas;
        t=cos((2*m-1)*pi/2/Nmeas);
        t=(t+1)/2;
end


Y=A1.*sin(2*pi*t)+A2.*cos(2*pi*t)+eps;
SSres_obs=getSSres(Y,t);


[~,I]=sort(permMat,2);

zts=t;
x1=gpuArray(sin(2*pi*t));
x2=gpuArray(cos(2*pi*t));
x0=gpuArray(ones(1,size(Y,2)));
X= [x0' x1' x2'];

betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y(I))));
fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-Y(I)).^2,2));
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
