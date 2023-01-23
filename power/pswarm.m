close all
clear
param.useGPU=true;
param.NL=5;
param.NR=3;
param.per=2; % period used in regression model
param.Amp=3;

simtype='fast';
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e2; % SMALL RIGHT NOW
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

costfun=@(t) -min(simulatePWR(param,cumsum(abs(t))/sum(abs(t))));

min(simulatePWR(param,'uniform'))
%%
swarm_opts=optimoptions(@particleswarm,'Display','Iter','MaxIterations',10);
topt=particleswarm(costfun,8,[],[],swarm_opts);
clf
xline(cumsum(abs(topt))/sum(abs(topt)))