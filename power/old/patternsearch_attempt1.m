clear
clf
simtype='medium';
solver='pattern_search'; % options simulanneal or pswarm
checkUseGPU

param.NL=4;
param.NR=4;
param.freq_true=7.5; % freq used in regression model
param.Amp=1;
param.noise1=1;
param.noise2=.5; % the noise actually used in the simulation
Nacro=32;
param.Nacro=Nacro;
nodes='uniform';
hold on
ylim([0,1])

repar=@(t) cumsum(abs(t))/sum(abs(t));
switch solver
    case 'pswarm'
        swarm_opts=optimoptions(@particleswarm,'Display','Iter','MaxIterations',2,'UseParallel',true);
    case {'simulanneal', 'sa_no_repar'}
        swarm_opts=optimoptions(@simulannealbnd,'Display','Iter','MaxIterations',100);
    case 'pattern_search'
        swarm_opts=optimoptions(@patternsearch,'Display','Iter', ...
            'MeshTolerance',1e-3,'FunctionTolerance',1e-3,'MaxIterations',100,'UseParallel',true);
end

acrovec=linspace(0,2*pi,Nacro+1);
acrovec=acrovec(1:end-1);

[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5);

if strcmp(solver,'sa_no_repar')
    costfun=@(t) -min(simulatePWR(param,sort(t)));
else
    costfun=@(t) -min(simulatePWR(param,repar(t)));
end

% get uniform value for reference
p_unif=simulatePWR(param,'uniform');
    
% optimize and plot result
switch solver
    case 'pswarm'
        [topt,fval]=particleswarm(costfun,param.NL+param.NR,[],[],swarm_opts);
    case 'simulanneal'
        [topt,fval]=simulannealbnd(costfun,t_unif,[],[],swarm_opts);
    case 'sa_no_repar'
        [topt,fval]=simulannealbnd(costfun,t_unif,zeros(1,param.NL+param.NR),ones(1,param.NL+param.NR),swarm_opts);
    case 'pattern_search'
        [topt,fval]=patternsearch(costfun,t_unif,[],[],[],[],[],[],[],swarm_opts);
end

if strcmp(solver,'sa_no_repar')
    topt_repar=topt;
else
    topt_repar=repar(topt);
end

simtype='long';

checkUseGPU
p_unif=simulatePWR(param,topt_repar);
p_nu=simulatePWR(param,topt_repar);
nexttile(1)
plot(acrovec,p_nu,'--b')
hold on
plot(acrovec,p_unif,'-k')
ylim([0,1])
hold on
drawnow
