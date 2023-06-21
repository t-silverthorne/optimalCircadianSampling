addpath('utils_core/')
addpath('utils_cost_fun/')

parpool_size=6;max_iter=1;Nbatch=1;

p.permMethod       = 'naive_reuse_perms';% p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
p.Nmeas     = 8;
p.Nacro     = 32;
p.Nresidual = 1e3;
p.Nperm     = 2e2;
p.noise     = 1;
p.Nbatch    = Nbatch;
p.Amp       = 2;
p.freq      = 3.8;

eps_cstr=5e-3;
Aineq=eye(p.Nmeas-1,p.Nmeas); 

for ii=1:p.Nmeas-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(p.Nmeas-1,1);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);

scale = wrap_getCostFun(t_unif,p,1:5);

opts=optimoptions('fminimax','Display','iter',...
                      'UseParallel',false,'MaxIterations',5);


[xopt_pareto,fopt_pareto] = fminimax(@(t) wrap_getCostFun(t,p,1:5)./abs(scale), ...
    sort(rand(1,p.Nmeas)),Aineq,bineq,[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],opts)