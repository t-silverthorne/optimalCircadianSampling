function pwr=simulatePWR(param,nodes)
method='QR'; % options are QR (uses matlab backslash) or normal (uses pseudo-inverse)

NL=param.NL;
NR=param.NR;
Nperm=param.Nperm;
Nresidual=param.Nresidual;

if ~isfield(param,'noise2')
    param.noise2=1;
end

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

pwr=[];
acrovec=linspace(0,2*pi,Nacro);
for i=1:Nacro
acro=acrovec(i);

if param.useGPU
    permMat=rand(Nresidual,Nmeas,Nperm,'gpuArray');
    eps=randn(Nresidual,Nmeas,'gpuArray');
else
    permMat=rand(Nresidual,Nmeas,Nperm);
    eps=randn(Nresidual,Nmeas);
end
[~,I]=sort(permMat,2);
Y=Amp*cos(2*pi*t*freq_true-acro)+param.noise2*eps;

X=constructX(t,param); % construct linear model

switch method
    case 'QR'
        betas_obs=X\Y'; % observed error
    case 'normal'
        betas_obs=(X'*X)\(X'*Y'); % observed error
end

fits_obs=(X*betas_obs)';
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));

m=size(Y,1);n=size(Y,2);
offMat=repmat((0:m-1)',1,n)*n;
Yp=Y';
YI=pagetranspose(Yp(pagetranspose(I+offMat)));

switch method
    case 'QR'
        betas=pagemldivide(X,pagetranspose(YI));
    case 'normal'
        betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
end
fits =pagetranspose(pagemtimes(X,betas));
SSres=sqrt(sum((fits-YI).^2,2));
pwr(end+1)=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/Nresidual;
end
end