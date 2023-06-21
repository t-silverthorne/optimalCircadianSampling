function [acrovec,pwr] = simulatePWR_slow(param,nodes)
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

pwr=[];
acrovec=linspace(0,2*pi,Nacro);
X=constructX(t,param); % construct linear model

for ii=1:Nacro
    acro=acrovec(ii);
    resstat=[];
    for jj=1:Nresidual
        eps=randn(1,Nmeas);
        Y=Amp*cos(2*pi*t*freq_true-acro)+param.noise*eps;
        
        betas_obs=(X'*X)\(X'*Y'); % observed error
        
        fits_obs=(X*betas_obs)';
        SSres_obs=sqrt(sum((fits_obs-Y).^2,2));
        permstat=[];
        for kk=1:Nperm
            Ypermuted=Y(randperm(Nmeas));
            
            betas=pagemldivide(X'*X,pagemtimes(X',pagetranspose(Ypermuted)));
            
            fits =pagetranspose(pagemtimes(X,betas));
            SSres=sqrt(sum((fits-Ypermuted).^2,2));
            permstat(end+1)=SSres>SSres_obs;
        end
        resstat(end+1)=sum(permstat)/Nperm; % on average how often was the true fit better
    end
    pwr(end+1)=sum(resstat>.95)/Nresidual;
end
end

