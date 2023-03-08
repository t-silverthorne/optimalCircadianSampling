function [acrovec,pwr]=simulatePWR_fully_vectorized(param,nodes)
NL=param.NL;
NR=param.NR;
Nperm=param.Nperm;
Nresidual=param.Nresidual;
Nacro=param.Nacro;
Amp=param.Amp;
Nmeas=NL+NR;
freq_true=param.freq_true; % do not confuse with freq_est

if isnumeric(nodes)
    t=nodes;
    Nmeas=length(t);
else
    switch nodes
        case 'uniform'
            [t,~]=getSamplingSchedules(NL,NR,0,0.5);
        case 'non-uniform'
            [~,t]=getSamplingSchedules(NL,NR,0,0.3);
    end

end

acrovec=linspace(0,2*pi,Nacro);
acromat=reshape(acrovec,1,1,1,Nacro);

if param.useGPU
    eps=randn(Nresidual,Nmeas,1,Nacro,'gpuArray');
else
    eps=randn(Nresidual,Nmeas,1,Nacro);
end
use_same_perms=true;

if use_same_perms
    if param.useGPU
        permMat=rand(Nresidual,Nmeas,Nperm,'gpuArray');
    else
        permMat=rand(Nresidual,Nmeas,Nperm);
    end
    [~,I]=sort(permMat,2);
    I=repmat(I,1,1,1,Nacro);
else
    if param.useGPU
        permMat=rand(Nresidual,Nmeas,Nperm,Nacro,'gpuArray');
    else
        permMat=rand(Nresidual,Nmeas,Nperm,Nacro);
    end

    [~,I]=sort(permMat,2);
end

clear permMat
Y=Amp*cos(2*pi*t*freq_true-acromat)+param.noise*eps;
clear eps

X=constructX(t,param); % construct linear model
if ~param.useGPU
    X=gather(X);
end

betas_obs=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error

fits_obs=pagetranspose(pagemtimes(X,betas_obs));
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));
clear fits_obs

m=size(Y,1);n=size(Y,2);r=size(Y,3);s=size(Y,4);
inds=0:m*n*r:m*n*r*s;
inds=inds(1:end-1);
inds=reshape(inds,1,1,1,length(inds));
offMat=repmat((0:m-1)',1,n)*n;
Yp=pagetranspose(Y);
YI=pagetranspose(Yp(pagetranspose(I+offMat+inds)));
clear I offMat inds
%max(max(max(max(abs(sort(YI,2)-sort(Y,2))))))

betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));

fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nresidual;
pwr=reshape(pwr,1,Nacro);
end