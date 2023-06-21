function [acrovec,pwr]=simulatePWR_matperm(param,nodes)
method='normal'; % options are QR (uses matlab backslash) or normal (uses pseudo-inverse)

NL=param.NL;
NR=param.NR;

Nresidual=param.Nresidual;

Nacro=param.Nacro;
Amp=param.Amp;
Nmeas=NL+NR;

Nperm=min(param.Nperm,factorial(Nmeas)); % sampling without replacement
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

pwr=[];
acrovec=linspace(0,2*pi,Nacro);
for i=1:Nacro
acro=acrovec(i);

if param.useGPU
    eps=randn(Nresidual,Nmeas,'gpuArray');
else
    eps=randn(Nresidual,Nmeas);
end


Y=Amp*cos(2*pi*t*freq_true-acro)+param.noise*eps;
X=constructX(t,param); % construct linear model

betas_obs=(X'*X)\(X'*Y'); % observed error
fits_obs=(X*betas_obs)';
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));


pall=perms(1:Nmeas);
Pmat=NaN(Nmeas,Nmeas,factorial(Nmeas));
I=eye(Nmeas);
for ii=1:factorial(Nmeas)
    Pmat(:,:,ii)= I(:,pall(ii,:)); % might be possible to vectorize
end

inds=randsample(1:factorial(Nmeas),Nperm,false);
YI=pagemtimes(Y,Pmat(:,:,inds));
betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));

fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr(end+1)=sum(sum(SSres>SSres_obs,3)/factorial(Nmeas)>.95)/Nresidual;
end
end
