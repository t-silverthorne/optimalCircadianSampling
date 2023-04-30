% popu_size also needs to be set at script call

testing=false;
addpath('utils_core/')
addpath('utils_cost_fun/')
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')
p.permMethod       = 'FY_double_for_loop'; p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'

if testing
    popu_size   = 3;
    p.Nmeas     = 8;
    p.Nacro     = 4;
    p.Nresidual = 1e1;
    p.Nperm     = 1e1;
    p.noise     = 1;
    p.Nbatch    = 1;
    p.Amp       = 2;
    p.freq      = 3.8;
else
    % needs parpool_size, max_iter, popu_size
    parpool(parpool_size); 
    p.Nmeas     = 8;
    p.Nacro     = 32;
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.noise     = 1;
    p.Nbatch    = 40;
    p.Amp       = 2;
    p.freq      = 3.8;
end

eps_cstr=5e-3;
Aineq=eye(p.Nmeas-1,p.Nmeas); 

for ii=1:p.Nmeas-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(p.Nmeas-1,1);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);

if testing
    opts=optimoptions(@patternsearch,'Display','iter','maxiter',2);
else
    opts=optimoptions(@patternsearch,'Display','iter','useparallel',true,'maxiter',max_iter);
end
x = patternsearch(@(t) -costfun_wrap_getPowerBatch(t,p),t_unif,Aineq,bineq,[],[],zeros(1,8),ones(1,8),[],opts) ;

nexttile(1)
[pwr,~,~]=getPowerBatch(t_unif,p);
pwr_unif=mean(pwr,1);
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1);
plot(acrovec,pwr_unif)
hold on

[pwr,~,~]=getPowerBatch(x,p);
pwr_opt=mean(pwr,1);
plot(acrovec,pwr_opt)
legend({'uniform','optimal'},'Interpreter','latex','location','southoutside','NumColumns',2)
%% pareto plot
% just worry about power and acro variance
if testing
    opts_pareto=optimoptions('paretosearch','Display','iter',...
                      'UseParallel',false, 'InitialPoints',repmat(t_unif,popu_size,1),...
                      'MeshTolerance',1e-2,'MaxIterations',3);
else
    opts_pareto=optimoptions('paretosearch','Display','iter',...
                      'UseParallel',true, 'InitialPoints',repmat(t_unif,popu_size,1),...
                      'MeshTolerance',1e-2,'MaxIterations',20);
end
nexttile(2)
Npar=5;
for ii=1:Npar
    [xopt_pareto,fopt_pareto] = paretosearch(@(t) wrap_getCostFun(t,p,[4,5]),p.Nmeas,Aineq,bineq,[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],opts_pareto);
    for ii=1:size(fopt_pareto,1)
        plot(fopt_pareto(ii,1),fopt_pareto(ii,2),'.k')
    end
    hold on
end
savefig('figs/fig2_bottom.fig')


function pmin=costfun_wrap_getPowerBatch(t,p)
[pwr,~,~]=getPowerBatch(t,p);
pmin=min(mean(pwr,1));
end
