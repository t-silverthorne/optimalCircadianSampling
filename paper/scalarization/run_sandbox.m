active_inds = 1:5;

nouter = 60;                    % number of outer optimzation loops
ninner = 1000;                    % number of inner optimization loops
parpool_size=30;
Nbatch=8;
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
p.Amp       = 8;
p.freq      = 3.8;
p.positive_cost_fun = true;

%% run particle swarm
eps_cstr=1e-2;


pmethod=2; % how should the design be parameterized

switch pmethod
    case 1
        xmaster=NaN(nouter,p.Nmeas-1);
    case 2
        xmaster=NaN(nouter,p.Nmeas);
end

parfor ii=1:nouter
    opts = optimoptions(@simulannealbnd,'Display','iter', ...
            'MaxIterations',ninner,'DisplayInterval',1,'ReannealInterval',50);
    
    switch pmethod
    	case 1
	    x0=rand(1,p.Nmeas-1);
	    xopt = simulannealbnd(@(dt) hvol_scalar_costfun(dt,p,active_inds,lambdaMaster(ii,:)), ...
				       x0,eps_cstr*ones(1,p.Nmeas-1),ones(1,p.Nmeas-1),opts);
	case 2
	    x0=rand(1,p.Nmeas);
	    xopt = simulannealbnd(@(t) hvol_scalar_costfun_v2(t,p,active_inds,lambdaMaster(ii,:)), ...
				       x0,eps_cstr*ones(1,p.Nmeas),ones(1,p.Nmeas),opts);
    end
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

function J=hvol_scalar_costfun_v2(t,p,active_inds,lambdaloc)
t=sort(t);
J=max(wrap_getCostFun(t,p,active_inds)./lambdaloc);
end
