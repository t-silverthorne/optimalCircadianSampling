% test of wrapper for cost function
p.Nmeas     = 8;
p.Nacro     = 16;
p.Nresidual = 2e2;
p.Nperm     = 2e2;
p.freq      = 3.9;
p.Amp       = 1.5;
p.noise     = 1.5;
p.Nbatch    = 1;
addpath('../utils_core')

clf
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[I3,I4]=constructUtilMats(p);
N=3*5;
sc=1e-1;    
wrap_getCostFun(t_unif,p,I3,I4)

% plot power only
for ii=1:N
    if ii<=N/3
        [pwr,~,~]=getPowerBatch(t_unif,p,I3,I4);
        plot(pwr,'-k')
    elseif N/3<ii && ii<=2*N/3
        [pwr,~,~]=getPowerBatch(rand(1,p.Nmeas),p,I3,I4);
        plot(pwr,'-r')    
    else
        [pwr,~,~]=getPowerBatch(t_unif+sc*rand(1,p.Nmeas),p,I3,I4);
        plot(pwr,'-b')    
    end
    hold on
    ylim([0,1])
    drawnow
end
% power histogram

% 
% Junif_mat=NaN(N,5);
% Jrand_mat=NaN(N,5);
% Jjit_mat =NaN(N,5);
% 
% parfor ii=1:N
%     Junif_mat(ii,:)=wrap_getCostFun(t_unif,p,I3,I4);
%     Jrand_mat(ii,:)=wrap_getCostFun(rand(1,p.Nmeas),p,I3,I4);
%     Jjit_mat(ii,:)=wrap_getCostFun(t_unif+sc*rand(1,p.Nmeas),p,I3,I4);
% end
% clf
% tiledlayout(5,1)
% nbins=max(10,floor(sqrt(N)));
% 
% for ii=1:5
%     nexttile(ii)
%     histogram(Jjit_mat(:,ii),nbins,'Normalization','probability','FaceColor','blue','EdgeColor','none')
%     hold on
%     histogram(Junif_mat(:,ii),nbins,'Normalization','probability','FaceColor','black','EdgeColor','none')
%     histogram(Jrand_mat(:,ii),nbins,'Normalization','probability','FaceColor','red','EdgeColor','none')
% end
% 
% 
% 
% %%
% 
% [t,~]=getSamplingSchedules(p.Nmeas,0,0,0);
% [I3,I4]=constructUtilMats(p);
% tic
% wrap_getCostFun(t,p,I3,I4)
% toc
% %%
% tic
% getPowerBatch(t,p,X,I3,I4);
% toc