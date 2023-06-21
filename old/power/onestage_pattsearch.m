close all
clear
simtype='fast';          % how many Monte Carlo samples to include
solver='pattern_search'; % options simulanneal or pswarm
addpath('utils/')
checkUseGPU              % uses simType to construct param
param.NL=4;              % samples in left interval
param.NR=4;              % samples in right interval
Nmeastot=param.NL+param.NR;
param.freq_true=8.2;     % freq used in regression model
param.Amp=3;             % signal amplitude
param.noise=1;         % the noise actually used in the simulation


Nacro=32;                % set up acrophase grid
param.Nacro=Nacro;
acrovec=linspace(0,2*pi,Nacro+1);
acrovec=acrovec(1:end-1);

% options for patternsearch
swarm_opts=optimoptions("patternsearch",'Algorithm',"nups",... 
            'Display','Iter', ...
            'MeshTolerance',1e-3,'MaxIterations',1000,'UseParallel',true, ...
            'PlotFcn',@show_current_schedule);

%'MeshContractionFactor',0.9,'MeshExpansionFactor',4.0, ...

Aineq=eye(Nmeastot-1,Nmeastot); % inequality constraints
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=ones(Nmeastot-1,1);

[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

% function to be optimized
costfun=@(t) -min(wrapsimulatePWR(param,t));
   
[topt,fval]=patternsearch(costfun,t_unif,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],swarm_opts);
% make a nice plot
close all

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

simtype='medium';
checkUseGPU
[acrovec,p_unif]=simulatePWR(param,t_unif);
[~,p_nu]=simulatePWR(param,topt);
nexttile(1)
hold on
ylim([0,1])

plot(acrovec,p_unif,'-k')
hold on
yline(min(p_unif),'--k')
plot(acrovec,p_nu,'-b')
yline(min(p_nu),'--b')

ylim([0,1])
hold on
drawnow

ylabel('$\beta(\phi;\xi)$','interpreter','latex')
xlabel('$\phi$','interpreter','latex')
xlim([0,2*pi])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
set(findall(gcf,'-property','FontSize'),'FontSize',10)

% legend
legend({'$\beta(\phi;\xi_{\mathrm{unif}})$', ...
    '$\mathrm{min}_\phi(\beta(\phi;\xi_{\mathrm{unif}}))$', ...
    '$\beta(\phi;\xi_{\mathrm{opt}})$', ...
    '$\mathrm{min}_\phi(\beta(\phi;\xi_{\mathrm{opt}}))$', ...
    },'Location','north','NumColumns',2,'FontSize',9)
plot_filename='figs/fig1_opt'
ht=3; % height
wd=4.5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too


function stop = show_current_schedule(x, optimValues, state)
stop = false;
clf
xlim([0,1])
xline(x.x)
hold on

drawnow
end
