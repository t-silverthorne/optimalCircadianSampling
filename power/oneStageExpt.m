clear
clf
tic
simtype='rough';
solver='sa_no_repar'; % options simulanneal or pswarm
checkUseGPU
tic
param.NL=4;
param.NR=4;
param.freq_true=7.5; % freq used in regression model
param.Amp=1;
Nacro=32;
param.Nacro=Nacro;
param.noise1=1;
param.noise2=1; % the noise actually used in the simulation
nodes='uniform';
hold on
ylim([0,1])

noisevals=logspace(log10(1),log10(0.4),3);
%noisevals=1/1.5;
%repar=@(t) 0.5*(1+tanh(sort(t)));
repar=@(t) cumsum(abs(t))/sum(abs(t));
switch solver
    case 'pswarm'
        swarm_opts_psw=optimoptions(@particleswarm,'Display','Iter','MaxIterations',2,'UseParallel',true);
    case {'simulanneal', 'sa_no_repar'}
        swarm_opts_sab=optimoptions(@simulannealbnd,'Display','Iter','MaxIterations',100);
end

acrovec=linspace(0,2*pi,Nacro+1);
acrovec=acrovec(1:end-1);
Nfreq=1;
tiledlayout(1,Nfreq,'TileSpacing','tight','Padding','tight')
[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5);

cvals=linspace(0.8,0,length(noisevals))'.*[1 1 1];
fprintf('running\n')
for ii=1:length(noisevals)
    fprintf('%d\n',ii)
    param.noise2=noisevals(ii);
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
            [topt,fval]=particleswarm(costfun,param.NL+param.NR,[],[],swarm_opts_psw);
        case 'simulanneal'
            %[topt,fval]=simulannealbnd(costfun,atanh(2*t_unif-1),[],[],swarm_opts_sab);
            [topt,fval]=simulannealbnd(costfun,t_unif,[],[],swarm_opts_sab);
        case 'sa_no_repar'
            [topt,fval]=simulannealbnd(costfun,t_unif,zeros(1,param.NL+param.NR),ones(1,param.NL+param.NR),swarm_opts_sab);
    end
    topt=sort(topt);
    if strcmp(solver,'sa_no_repar')
        topt_repar=topt;
    else
        topt_repar=repar(topt);
    end

    p_nu=simulatePWR(param,topt_repar);
    nexttile(1)
    plot(acrovec,p_nu,'Color',cvals(ii,:))
    hold on
    plot(acrovec,p_unif,'--k','Color',cvals(ii,:))
    ylim([0,1])
    hold on
%     nexttile(Nfreq+1)
%     ylim([-1,1])
%     plot(acrovec,p_nu-p_unif,'Color',cvals(ii,:))
    hold on
    drawnow
end
toc