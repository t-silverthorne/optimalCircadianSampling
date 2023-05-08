testing=false;
parpool(parpool_size)
addpath('utils_core/')
addpath('utils_cost_fun/')
tiledlayout(2,5,'TileSpacing','tight','Padding','tight')

p.permMethod='fast'
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

eps_cstr=5e-3;
Aineq=eye(p.Nmeas-1,p.Nmeas); 

for ii=1:p.Nmeas-1
    Aineq(ii,ii+1)=-1;
end
bineq=-eps_cstr*ones(p.Nmeas-1,1);
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);

if testing
    opts_pareto=optimoptions('paretosearch','Display','iter',...
                      'UseParallel',false, 'InitialPoints',repmat(t_unif,popu_size,1),...
                      'MeshTolerance',1e-2,'MaxIterations',max_iter);
else
    opts_pareto=optimoptions('paretosearch','Display','iter',...
                      'UseParallel',true, 'InitialPoints',repmat(t_unif,popu_size,1),...
                      'MeshTolerance',1e-2,'MaxIterations',max_iter);
end



% phivec=linspace(0,2*pi,p.Nacro+1);
% phivec=phivec(1:end-1);
fopt_pareto_master=[]
%while size(fopt_pareto_master,1)<5e2
while size(fopt_pareto_master,1)<num_pareto_points
    [xopt_pareto,fopt_pareto] = paretosearch(@(t) wrap_getCostFun(t,p,I3,I4,[1:5]),p.Nmeas,Aineq,bineq,[],[],zeros(1,p.Nmeas),ones(1,p.Nmeas),[],opts_pareto);
    fopt_pareto_master=[fopt_pareto_master; fopt_pareto];
end


for ii=1:5
    nexttile(ii)
    histogram(fopt_pareto_master(:,ii))
end

for ii=1:4
    nexttile(5+ii)
    for jj=1:size(fopt_pareto_master,1)
        plot(fopt_pareto_master(jj,ii),fopt_pareto_master(jj,5),'.k')
        hold on
    end
end
savefig('figs/fig3.fig')
