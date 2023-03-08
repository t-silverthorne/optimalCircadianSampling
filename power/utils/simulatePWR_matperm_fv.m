function [acrovec,pwr] = simulatePWR_matperm_fv(param,nodes)
NL=param.NL;
NR=param.NR;

Nresidual=param.Nresidual;
Nacro=param.Nacro;
Amp=param.Amp;
Nmeas=NL+NR;

Nperm=min(param.Nperm,factorial(Nmeas));

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
Y=Amp*cos(2*pi*t*freq_true-acromat)+param.noise*eps;
clear eps

X=constructX(t,param); % construct linear model
if ~param.useGPU
    X=gather(X);
end

betas_obs=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error

fits_obs=pagetranspose(pagemtimes(X,betas_obs));
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));
clear betas_obs fits_obs

if ~isfield(param,'Pmat')
    pall=perms(1:Nmeas);
    Pmat=NaN(Nmeas,Nmeas,factorial(Nmeas));
    I=eye(Nmeas);
    for ii=1:factorial(Nmeas)
        Pmat(:,:,ii)= I(:,pall(ii,:)); % might be possible to vectorize
    end
    param.Pmat=Pmat;
    clear Pmat
end
inds=randsample(1:factorial(Nmeas),Nperm,false);
YI=pagemtimes(Y,param.Pmat(:,:,inds));

betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));

fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nresidual;
pwr=reshape(pwr,1,Nacro);
end
