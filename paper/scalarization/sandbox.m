% generate weights
active_inds = 1:5;

nouter = 10; % number of outer optimzation loops
ninner = 10; % number of inner optimization loops
parpool_size=6;Nbatch=5;popu_size=24;

numOpt = length(active_inds);    % dimension of multi-objective optimization problem
% construct weight matrix
addpath('../utils_core/')
addpath('../utils_cost_fun/')

lambdaMaster=[];
while size(lambdaMaster,1)<nouter
    lambdaMat=randn(nouter,numOpt);
    lambdaMat=lambdaMat./sqrt(sum(lambdaMat.^2,2));
    lambdaMaster = [lambdaMaster; lambdaMat(sum(lambdaMat>0,2)==numOpt,:)];
end
lambdaMaster
% run particle swarm

p.permMethod= 'naive_reuse_perms';% p.permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'
p.Nmeas     = 8;
p.Nacro     = 32;
p.Nresidual = 1e3;
p.Nperm     = 2e2;
p.noise     = 1;
p.Nbatch    = Nbatch;
p.Amp       = 2;
p.freq      = 3.8;
p.positive_cost_fun = true;

eps_cstr=1e-2;

xmaster=NaN(nouter,p.Nmeas-1);

%%
for ii=1:nouter
    % construct swarm matrix
    A=rand(popu_size,p.Nmeas-1);
    SM=A./(.1*rand(popu_size,1)+sum(A,2));
    
    % optimization options
    opts = optimoptions(@particleswarm,'Display','iter','SwarmSize', popu_size, ...
            'UseParallel',true,'MaxIterations',ninner,'InitialSwarmMatrix',SM);
    xopt = particleswarm(@(dt) hvol_scalar_costfun(dt,p,active_inds,lambdaMaster(ii,:)),p.Nmeas-1,eps_cstr*ones(1,p.Nmeas-1),ones(1,p.Nmeas-1)/p.Nmeas,opts);
    
    xmaster(ii,:)=xopt;
end

save('scalar_testrun.mat')
%%%
%disp("finished1")
%clf
%for ii=1:nouter
%    t=[0 cumsum(xmaster(ii,:))];
%    J=wrap_getCostFun(t,p,active_inds)
%    plot(J(4),J(5),'ok')
%    hold on
%end
%disp("finished2")














function J=hvol_scalar_costfun(dt,p,active_inds,lambdaloc)
t=[0 cumsum(dt)];
J=max(wrap_getCostFun(t,p,active_inds)./lambdaloc);
end
