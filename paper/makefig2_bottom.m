parpool(parpool_size); 
% popu_size also needs to be set at script call

testing=false;
addpath('utils_core/')
addpath('utils_cost_fun/')
tiledlayout(2,1,'TileSpacing','tight','Padding','tight')
p.permMethod       = 'FY_double_for_loop'; permActionMethod = 'index'; % options index or matrix for 'naive_make_perms_first'

if testing
    p.Nmeas     = 8;
    p.Nacro     = 4;
    p.Nresidual = 1e1;
    p.Nperm     = 1e1;
    p.noise     = 1;
    p.Nbatch    = 1;
    p.Amp       = 2;
    p.freq      =3.8;
else
    p.Nmeas     = 8;
    p.Nacro     = 32;
    p.Nresidual = 1e3;
    p.Nperm     = 1e2;
    p.noise     = 1;
    p.Nbatch    = 20;
    p.Amp       = 2;
    p.freq      =3.8;
end
[I3,I4]=constructUtilMats(p);

% p.permMethod='fast';
% %% particle swarm
% opts = optimoptions(@particleswarm, ... %'HybridFcn',@fminsearch,...
%                                        'Display','iter', ...    
%                                        'SwarmSize',24, ...
%                                        'MaxIter',10, ... % TODO make realistic
%                                        'FunctionTolerance',1e-3, ... % might be too high
%                                        'UseParallel',true);%,'PlotFcn',@pswmyfun); % TODO make true
% 
% [t_opt,~]=particleswarm( @(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),8,zeros(1,8),ones(1,8),opts)
% %% fmincon

% [pwr,~,~]=getPowerBatch(t_unif,p,I3,I4)
% min(mean(pwr,1))
% 
% opts=optimoptions(@fmincon,'Display','iter','UseParallel',true,'MaxIter',10);
% x = fmincon(@(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),t_unif,Aineq,bineq,[],[],zeros(1,8),ones(1,8),[],opts);
% 
% %% fminsearch
% opts=optimset('Display','iter','MaxIter',10);
% x=fminsearch(@(t) -fminsearch_wrap_getPowerBatch(t,p,I3,I4),log10(t_unif),opts);

% patternsearch
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
    opts=optimoptions(@patternsearch,'Display','iter','useparallel',true,'maxiter',50);
end
x = patternsearch(@(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),t_unif,Aineq,bineq,[],[],zeros(1,8),ones(1,8),[],opts) ;

nexttile(1)
[pwr,~,~]=getPowerBatch(t_unif,p,I3,I4);
pwr_unif=mean(pwr,1);
acrovec=linspace(0,2*pi,p.Nacro+1);
acrovec=acrovec(1:end-1);
plot(acrovec,pwr_unif)
hold on

[pwr,~,~]=getPowerBatch(x,p,I3,I4);
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
    [xopt_pareto,fopt_pareto] = paretosearch(@(t) wrap_getCostFun(t,p,I3,I4,[4,5]),p.Nmeas,Aineq,bineq,[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],opts_pareto);
    for ii=1:size(fopt_pareto,1)
        plot(fopt_pareto(ii,1),fopt_pareto(ii,2),'.k')
    end
    hold on
end
savefig('figs/fig2_bottom.fig')

% nexttile(tileind(ny1+ny2+1,1),[ny3,nx1])
% plot(1,1)
% nexttile(tileind(ny1+ny2+ny3+1,1),[ny3,nx1])
% plot(1,1)
% %% sweeps
% nexttile(tileind(ny1+ny2+1,nx1+1),[ny3,nx2])
% plot(1,1)
% nexttile(tileind(ny1+ny2+ny3+1,nx1+1),[ny3,nx2])
% plot(1,1)
% 

%% global search
% problem = createOptimProblem('fmincon',...
%     'objective',@(t) -costfun_wrap_getPowerBatch(t,p,I3,I4),...
%     'x0',t_unif,'options',...
%     optimoptions(@fmincon,'Display','iter','UseParallel',false,'MaxIter',10),);


function pmin=arrayfun_wrap_getPowerBatch(Amp,freq,t,p,I3,I4)
p.Amp=p.noise*Amp;
p.freq=freq;
[pwr,~,~]=getPowerBatch(t,p,I3,I4);
pmin=min(mean(pwr,1));
end

function pmin=costfun_wrap_getPowerBatch(t,p,I3,I4)
[pwr,~,~]=getPowerBatch(t,p,I3,I4);
pmin=min(mean(pwr,1));
end

function pmin=fminsearch_wrap_getPowerBatch(t,p,I3,I4)
t=sort(0.5*(1+tanh(10.^t)));
pmin=costfun_wrap_getPowerBatch(t,p,I3,I4);
end
