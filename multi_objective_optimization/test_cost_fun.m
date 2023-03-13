clear
gpuDevice(1)
addpath('utils/')
simtype='medium';
checkUseGPU
param.useGPU=true;
param.useParallel=false;
param.NL=4;
param.NR=4;
param.Nmeas=param.NL+param.NR;
param.freq_true=7.6;
param.Amp=1;
param.noise=.5;
param.Nacro=16;
nodes='uniform';
Nmeastot=param.NL+param.NR;
[t_unif,~]=getSamplingSchedules(param.NL,param.NR,0,0.5); % initial guess for sampling

param.Pmat=construct_all_perms(param.Nmeas);
%%
% multi objective func to be minimized

tic
costfun('uniform')
toc

%%
eps_cstr=1e-2;
Aineq=eye(Nmeastot-1,Nmeastot); % inequality constraints
for ii=1:Nmeastot-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(Nmeastot-1,1);

%% Run constrained optimization


% paretosearch (pattern search for multiobjective
tic
opts=optimoptions('paretosearch','Display','iter',...
                  'UseParallel',param.useParallel, ...
                  'MeshTolerance',1e-2,'MaxIterations',5);
[xopt,~] = paretosearch(costfun,Nmeastot,Aineq,bineq,[],[],zeros(1,Nmeastot),ones(1,Nmeastot),[],opts);
toc
test_result(param,xopt)
%% gamultiobj (genetic algorithm for multiobjective)
param.useParallel=true;
param.useGPU=~param.useParallel;
Nmeastot=param.NL+param.NR;
costfun = @(t) costfun_power_bias_var(param,t);
tic
param.useParallel=false;
opts=optimoptions('gamultiobj','Display','iter','InitialPopulationMatrix',t_unif,'UseParallel',param.useParallel,...
                      'PopulationSize',50,'MaxGenerations',5);
[xopt,~] = gamultiobj(costfun,Nmeastot,Aineq,bineq,[],[],zeros(Nmeastot,1),ones(Nmeastot,1),[],opts);
toc

test_result(param,xopt)
% Run reparameterized optimization


%% Penalty methods

function test_result(param,xopt)
% Compare uniform solution, random solution and all pareto front solutions
Nmeastot=param.NL+param.NR;
clf

tiledlayout(5,3,'TileIndexing','columnmajor')
[acrovec,pwr,est]=wrap_simulatePWR_matperm_fv(param,'uniform');
plot_multi_obj(param,acrovec,pwr,est,'-k',1)
drawnow

for ii=1:size(xopt,1)
    [acrovec,pwr,est]=wrap_simulatePWR_matperm_fv(param,sort(rand(1,Nmeastot)));
    plot_multi_obj(param,acrovec,pwr,est,'-b',2)
end
drawnow

for ii=1:size(xopt,1)
    [acrovec,pwr,est]=wrap_simulatePWR_matperm_fv(param,xopt(ii,:));
    plot_multi_obj(param,acrovec,pwr,est,'-r',3)
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
    for jj=1:3
        h=nexttile(5*(jj-1)+ii);
        h.YScale='log'
        h.YLim=[y1min y2max];
    end
end
end


function plot_multi_obj(param,acrovec,pwr,est,lt,ind)

% power
nexttile(5*(ind-1)+1)
hold on
plot(acrovec,pwr,lt)

% amp bias
nexttile(5*(ind-1)+2)
hold on
semilogy(acrovec,abs(param.Amp-est.amp_mu),lt)

% amp variance
nexttile(5*(ind-1)+3)
hold on
semilogy(acrovec,abs(est.amp_st),lt)

% phase error
nexttile(5*(ind-1)+4)
hold on
semilogy(acrovec,abs(exp(1j*acrovec)-exp(1j*est.phi_mu)),lt)


% phase variance
nexttile(5*(ind-1)+5)
hold on
semilogy(acrovec,est.phi_cvar,lt)
end