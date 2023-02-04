clear
clf
simtype='fast';
solver='pattern_search'; % options simulanneal or pswarm
checkUseGPU

param.NL=4;
param.NR=4;
Nmeastot=param.NL+param.NR;
param.freq_true=8.5; % freq used in regression model
param.Amp=1;
param.noise1=0.7;
param.noise2=.5; % the noise actually used in the simulation
Nacro=32;
param.Nacro=Nacro;

swarm_opts=optimoptions("patternsearch",'Algorithm',"nups-mads",'Display','Iter', ...
            'MeshTolerance',5e-3,'MaxIterations',100,'UseParallel',true);
acrovec=linspace(0,2*pi,Nacro+1);
acrovec=acrovec(1:end-1);
[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5);

costfun=@(t) -min(simulatePWR(param,t));
   
Aineq=eye(Nmeastot-1,Nmeastot);
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=zeros(Nmeastot-1,1);


[topt,fval]=patternsearch(costfun,t_unif,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],swarm_opts);
%%
clf

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


simtype='long';
checkUseGPU
p_unif=simulatePWR(param,t_unif);
p_nu=simulatePWR(param,topt);
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
%%

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