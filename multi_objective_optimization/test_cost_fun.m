clear
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% things to tweak on CC
plot_filename = 'test_cost_fun';
popu_size     = 10;
max_gen       = 2;

addpath('utils/')
simtype='medium';
checkUseGPU
param.useGPU=false;
param.useParallel=false;

% waveform
param.NL=4;
param.NR=4;
param.Nmeas=param.NL+param.NR;
param.freq_true=7.4;
param.Amp=2.5;
param.noise=1;

% sample size
param.Nperm=1e2;
param.Nresidual=1e2;
param.Nacro=16; % num. fourier samples
param.method='4tensor';
param.perm_method='fy'; % options fy or randperm

% initial guess
nodes='uniform';
Nmeastot=param.NL+param.NR;
[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

% tic
% [~,p,~]=simulatePWR_rank4(param,'uniform')
% toc


% inequality constraints
eps_cstr=5e-3;
Aineq=eye(Nmeastot-1,Nmeastot); 
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(Nmeastot-1,1);

% optimization settings
param.useParallel=true;
param.useGPU=false;
Nmeastot=param.NL+param.NR;

costfun = @(t) costfun_power_bias_var(param,t); 

opts_GA=optimoptions('gamultiobj','Display','iter','UseParallel',param.useParallel,...
                      'PopulationSize',popu_size,'MaxGenerations',max_gen);
% run optimization
tic
[xopt_GA,fopt_GA] = gamultiobj(costfun,Nmeastot,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],opts_GA);
toc

tic
opts_pareto=optimoptions('paretosearch','Display','iter',...
                  'UseParallel',param.useParallel, 'InitialPoints',repmat(t_unif,popu_size,1),...
                  'MeshTolerance',1e-2,'MaxIterations',max_gen);
[xopt_pareto,fopt_pareto] = paretosearch(costfun,Nmeastot,Aineq,bineq,[],[],zeros(1,Nmeastot),ones(1,Nmeastot),[],opts_pareto);
toc

%% histogram plot
[acrovec,pwr,est]=simulatePWR_rank4(param,'uniform');

J_unif=costfun_power_bias_var(param,t_unif)

close all
nbin=max(10,sqrt(popu_size));
tiledlayout(3,5)
for ii=1:5
    nexttile(ii)
    histogram(fopt_GA(:,ii),nbin,Normalization="count")
    hold on
    histogram(fopt_pareto(:,ii),nbin,Normalization='count')
    xline(J_unif(ii))
end


%% plot uniform for comparison


%% plot jittered for comparison

savefig(gcf,strcat(plot_filename,'.fig'))% save matlab .fig too


%%
% test_result(param,xopt_GA,xopt_pareto)
% 
% 
% % clean graph
% nexttile(1)
% ylabel('Power')
% nexttile(2)
% ylabel('amp bias')
% nexttile(3)
% ylabel('amp std')
% nexttile(4)
% ylabel('phase bias')
% nexttile(5)
% ylabel('|order|')

function test_result(param,xopt1,xopt2)
% Compare uniform solution, random solution and all pareto front solutions
Nmeastot=param.NL+param.NR;
clf

tiledlayout(5,4,'TileIndexing','columnmajor')
[acrovec,pwr,est]=simulatePWR_rank4(param,'uniform');
plot_multi_obj(param,acrovec,pwr,est,'-k',1)
drawnow


for ii=1:size(xopt1,1)
    [acrovec,pwr,est]=simulatePWR_rank4(param,sort(rand(1,Nmeastot)));
    plot_multi_obj(param,acrovec,pwr,est,'-b',2)
end
drawnow


for ii=1:size(xopt1,1)
    [acrovec,pwr,est]=simulatePWR_rank4(param,xopt1(ii,:));
    plot_multi_obj(param,acrovec,pwr,est,'-r',3)
end
drawnow


for ii=1:size(xopt2,1)
    [acrovec,pwr,est]=simulatePWR_rank4(param,xopt2(ii,:));
    plot_multi_obj(param,acrovec,pwr,est,'-r',4)
end
hold off

% clean up axes
for ii=1:5
    y1min=Inf;
    y2max=-Inf;
    for jj=1:3
        h=nexttile(5*(jj-1)+ii);
        if h.YLim(1)<y1min
            y1min=h.YLim(1);
        end
        if h.YLim(2)>y2max
            y2max=h.YLim(2);
        end
        
    end
    if y1min==0
        y1min=1e-3;
    end
    for jj=1:4
        h=nexttile(5*(jj-1)+ii);
        h.YScale='log';
        h.YLim=[y1min y2max];
    end
end
end


function plot_multi_obj(param,acrovec,pwr,est,lt,ind)

% power
nexttile(5*(ind-1)+1)
hold on
plot(acrovec,pwr,lt)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})

% amp bias
nexttile(5*(ind-1)+2)
hold on
semilogy(acrovec,abs(param.Amp-est.amp_mu),lt)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})


% amp variance
nexttile(5*(ind-1)+3)
hold on
semilogy(acrovec,abs(est.amp_st),lt)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})

% phase error
nexttile(5*(ind-1)+4)
hold on
semilogy(acrovec,min(mod([acrovec-est.phi_mu; est.phi_mu-acrovec],2*pi),[],1),lt)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})


% phase variance
nexttile(5*(ind-1)+5)
hold on
semilogy(acrovec,est.phi_cvar,lt)
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'$$0$$','$$\frac{\pi}{2}$$','$$\pi$$','$$\frac{3\pi}{2}$$','$$2\pi$$'})

end