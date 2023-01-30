clear
clf
tic
simtype='medium';
checkUseGPU
tic
param.NL=4;
param.NR=4;
param.freq_true=1.5; % freq used in regression model
param.Amp=1;
param.noise1=1;
param.noise2=1; % the noise actually used in the simulation
nodes='uniform';
hold on
ylim([0,1])

noisevals=logspace(log10(2),log10(1/2),4);
costfun=@(t) -min(simulatePWR(param,cumsum(abs(t))/sum(abs(t))));
swarm_opts=optimoptions(@particleswarm,'Display','Iter','MaxIterations',5,'UseParallel',true);
acrovec=linspace(0,2*pi,param.Nacro);
Nfreq=1;
tiledlayout(2,Nfreq,'TileSpacing','tight','Padding','tight')

cvals=linspace(0.8,0,length(noisevals))'.*[1 1 1];
fprintf('running\n')
for ii=1:length(noisevals)
    fprintf('%d\n',ii)
    param.noise2=noisevals(ii);

    % get uniform value for reference
    p_unif=simulatePWR(param,'uniform');
    
    % optimize and plot result
    [topt,fval]=particleswarm(costfun,param.NL+param.NR,[],[],swarm_opts);
    p_nu=simulatePWR(param,topt);
    nexttile(1)
    plot(acrovec,p_nu,'Color',cvals(ii,:))
    ylim([0,1])
    hold on
    nexttile(Nfreq+1)
    ylim([-1,1])
    plot(acrovec,p_nu-p_unif,'Color',cvals(ii,:))
    hold on
    drawnow
end
toc