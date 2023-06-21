close all
clear

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

simtype='medium';          % how many Monte Carlo samples to include
addpath('utils/')
checkUseGPU              % uses simType to construct param
param.NL=4;              % samples in left interval
param.NR=4;              % samples in right interval
Nmeastot=param.NL+param.NR;
param.freq_true=8.1;     % freq used in regression model
param.Amp=2;             % signal amplitude
param.noise=1;          % the noise actually used in the simulation
param.Nacro=32;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling
[~,p_unif]=simulatePWR(param,t_unif);


t_unif_mod=t_unif(2:end);

costfun=@(t) -min(fminbndWrapSimulatePWR(param,t_unif_mod,0.5*(1+tanh(t))));
min(p_unif)

[t0,fval]=fminsearch(costfun,0);
fval
t_nu=[t_unif_mod 0.5*(1+tanh(t0))];
[acrovec,p_nu]=simulatePWR(param,t_nu);

plot(acrovec,p_unif,'-k')
hold on
plot(acrovec,p_nu,'-b')
ylim([0,1])
%%

ylabel('$\beta(\phi;\xi)$','interpreter','latex')
xlabel('$\phi$','interpreter','latex')
xlim([0,2*pi])
ylim([0,.2])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})
set(findall(gcf,'-property','FontSize'),'FontSize',10)

legend({'uniform','fminsearch'})
plot_filename='figs/only_one_point_moved'
ht=3; % height
wd=4.5; % width
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode','manual','PaperSize',[wd,ht],'PaperPosition',[0 0 wd ht])
print(gcf,plot_filename,'-dpng','-r600') % -r sets the resolution
savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too



function pwr = fminbndWrapSimulatePWR(param,nodes,t0)
nodes_full=[nodes t0];
pwr=wrapsimulatePWR(param,nodes_full);
end