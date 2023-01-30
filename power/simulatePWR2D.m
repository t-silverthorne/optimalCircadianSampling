clear
simtype='fast';
checkUseGPU
tic
param.NL=4;
param.NR=4;
param.freq_true=2.5; % freq used in regression model
param.acro_true='rand';
param.Amp=1;
param.noise1=1;
param.noise2=1; % the noise actually used in the simulation
nodes='uniform';

min(simulatePWR(param,nodes))
method='backslash'; % options are QR (uses matlab backslash) or normal (uses pseudo-inverse)

NL=param.NL;
NR=param.NR;
NMC=param.Nperm*param.Nresidual; % number of Monte Carlo samples

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
    permMat=rand(NMC,Nmeas,'gpuArray');
    eps=randn(NMC,Nmeas,'gpuArray');
else
    permMat=rand(NMC,Nmeas,Nperm);
    eps=randn(NMC,Nmeas);
end
[~,I]=sort(permMat,2);
Y=Amp*cos(2*pi*t*freq_true-acro)+param.noise2*eps;

X=constructX(t,param); % construct linear model

switch method
    case 'backslash'
        betas_obs=X\Y'; % observed error
    case 'normal'
        betas_obs=(X'*X)\(X'*Y'); % observed error
end

fits_obs=(X*betas_obs)';
SSres_obs=sqrt(sum((fits_obs-Y).^2,2));

m=size(Y,1);n=size(Y,2);
offMat=repmat((0:m-1)',1,n)*n;
Yp=Y';
YI=transpose(Yp(transpose(I+offMat)));

switch method
    case 'backslash'
        betas=X\transpose(YI);
    case 'normal'
        betas=(X'*X)\(X'*YI');
end
fits =(X*betas)';
SSres=sqrt(sum((fits-YI).^2,2));
pwr(end+1)=sum(sum(SSres>SSres_obs,3)/Nperm>.95)/NMC;
end