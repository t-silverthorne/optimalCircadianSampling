clear
p.Nmeas     = 8;
p.Nacro     = 16;
p.Nresidual = 5e2;
p.Nperm     = 2e2;
p.freq      = 7.8;
p.Amp       = 1;
p.noise     = .3;
p.Nbatch    = 1;

p.useParallel=true;
popu_size   =6*5;
max_gen     =5;

addpath('../utils_core')
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[X,I3,I4]=constructUtilMats(t_unif,p);
J_unif=wrap_getCostFun(t_unif,p,X,I3,I4)
%%
% inequality constraints
eps_cstr=5e-3;
Aineq=eye(p.Nmeas-1,p.Nmeas); 
for ii=1:p.Nmeas-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(p.Nmeas-1,1);

costfun = @(t) wrap_getCostFun(t,p,X,I3,I4); 

opts_GA=optimoptions('gamultiobj','Display','iter','UseParallel',p.useParallel,...
                      'PopulationSize',popu_size,'MaxGenerations',max_gen);
% run optimization
tic
[xopt_GA,fopt_GA] = gamultiobj(costfun,p.Nmeas,Aineq,bineq,[],[],zeros(p.Nmeas,1),ones(p.Nmeas,1),[],opts_GA);
toc

%%

tic
opts_pareto=optimoptions('paretosearch','Display','iter',...
                  'UseParallel',p.useParallel, 'InitialPoints',repmat(t_unif,popu_size,1),...
                  'MeshTolerance',1e-2,'MaxIterations',3);
[xopt_pareto,fopt_pareto] = paretosearch(costfun,p.Nmeas,Aineq,bineq,[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],opts_pareto);
toc

J_unif=wrap_getCostFun(t_unif,p,X,I3,I4);

close all
nbin=max(10,sqrt(popu_size));
tiledlayout(5,1)
for ii=1:5
    nexttile(ii)
    histogram(fopt_pareto(:,ii),nbin,Normalization="count")
    hold on
    xline(J_unif(ii))
end


