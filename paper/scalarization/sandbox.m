% generate weights
active_inds = 1:5;

nouter = 120/5;                    % number of outer optimzation loops
ninner = 100/5;                    % number of inner optimization loops
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

%% run particle swarm
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

eps_cstr=1e-2;

xmaster=NaN(nouter,p.Nmeas-1);

optmethod='simulannealbnd';
parfor ii=1:nouter

    switch optmethod
        case 'particleswarm'
            % construct swarm matrix
            A=rand(popu_size,p.Nmeas-1);
            SM=A./(.1*rand(popu_size,1)+sum(A,2));
            
            % optimization options
            opts = optimoptions(@particleswarm,'Display','iter','SwarmSize', popu_size, ...
                    'UseParallel',true,'MaxIterations',ninner, ...
                    'InitialSwarmMatrix',SM);
            xopt = particleswarm(@(dt) hvol_scalar_costfun(dt,p,active_inds,lambdaMaster(ii,:)), ...
                    p.Nmeas-1,eps_cstr*ones(1,p.Nmeas-1),ones(1,p.Nmeas-1)/p.Nmeas,opts);
        case 'simulannealbnd'
            opts = optimoptions(@simulannealbnd,'Display','iter', ...
                    'MaxIterations',ninner,'DisplayInterval',1,'ReannealInterval',15);
            x0=rand(1,p.Nmeas-1);
            xopt = simulannealbnd(@(dt) hvol_scalar_costfun(dt,p,active_inds,lambdaMaster(ii,:)), ...
                                       x0,eps_cstr*ones(1,p.Nmeas-1),ones(1,p.Nmeas-1)/p.Nmeas,opts);
    end
    xmaster(ii,:)=xopt;
end


save(fname)
%%
load('scalarized_test_simulanneal_larger_batch.mat')
close all
addpath('../utils_core/')
addpath('../utils_cost_fun/')
tiledlayout(3,2)

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

xlab_list={'amp bias','acro bias','amp variance', 'acro variance', '1-power'};


[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
J_unif = wrap_getCostFun(t_unif,p,active_inds);


disp("finished1")

for ii=1:size(xmaster,1)
    t=[0 cumsum(xmaster(ii,:))];
    J=wrap_getCostFun(t,p,active_inds);
    
    for jj=1:5
        nexttile(jj)    
        if jj==5
            plot(J(jj),J(1),'.k','MarkerSize',7);
            hold on
            if (ii==1 || ii == size(xmaster,1))
               plot(J_unif(jj),J_unif(1),'.r','MarkerSize',10);
            end
        else
            plot(J(jj),J(jj+1),'.k','MarkerSize',7);
            hold on
            if (ii==1 || ii == size(xmaster,1))
               plot(J_unif(jj),J_unif(jj+1),'.r','MarkerSize',20);
            end
        end
        
    end

    for jj=1:5
        nexttile(jj)
        xlabel(xlab_list{jj})
        if jj==5
            ylabel(xlab_list{1});
        else
            ylabel(xlab_list{jj+1});
        end
    end
    drawnow
   
end
disp("finished2")

function J=hvol_scalar_costfun(dt,p,active_inds,lambdaloc)
t=NaN(1,length(dt)+1);
t(1)=0;
for ii=2:length(t)
    t(ii) = t(ii-1) + dt(ii)*(1-t(ii-1));
end
J=max(wrap_getCostFun(t,p,active_inds)./lambdaloc);
end
