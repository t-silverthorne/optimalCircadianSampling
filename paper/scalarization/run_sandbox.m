active_inds = 1:5;

nouter = 120;                    % number of outer optimzation loops
ninner = 100;                    % number of inner optimization loops
parpool_size=6;
Nbatch=5;
popu_size=30;
fname='run_sandbox_out.mat';
parpool(parpool_size)
numOpt = length(active_inds);    % dimension of multi-objective optimization problem
% construct weight matrix
addpath('../utils_core/')
addpath('../utils_cost_fun/')


lambdaMaster=[];                 % Sample positive orthant of n sphere
while size(lambdaMaster,1)<nouter
    lambdaMat=randn(nouter,numOpt);
    lambdaMat=lambdaMat./sqrt(sum(lambdaMat.^2,2));
    lambdaMaster = [lambdaMaster; lambdaMat(sum(lambdaMat>0,2)==numOpt,:)];
end

p.permMethod= 'naive_reuse_perms';
p.Nmeas     = 8;
p.Nacro     = 16;
p.Nresidual = 1e3;
p.Nperm     = 1e2;
p.noise     = 1;
p.Nbatch    = Nbatch;
p.Amp       = 2;
p.freq      = 3.8;
p.positive_cost_fun = true;

%% run particle swarm

hvol_scalar_costfun(rand(1,p.Nmeas-1),p,active_inds,lambdaMaster(1,:))
%%
eps_cstr=1e-2;
xmaster=NaN(nouter,p.Nmeas-1);

parfor ii=1:nouter

    opts = optimoptions(@simulannealbnd,'Display','iter', ...
            'MaxIterations',ninner,'DisplayInterval',1,'ReannealInterval',15);
    x0=rand(1,p.Nmeas-1);
    xopt = simulannealbnd(@(dt) hvol_scalar_costfun(dt,p,active_inds,lambdaMaster(ii,:)), ...
                               x0,eps_cstr*ones(1,p.Nmeas-1),ones(1,p.Nmeas-1),opts);

    xmaster(ii,:)=xopt;
end

save(fname)

function J=hvol_scalar_costfun(dt,p,active_inds,lambdaloc)
t=NaN(1,length(dt)+1);
t(1)=0;
for ii=2:length(t)
    t(ii) = t(ii-1) + dt(ii-1)*(1-t(ii-1));
end
J=max(wrap_getCostFun(t,p,active_inds)./lambdaloc);
end
z