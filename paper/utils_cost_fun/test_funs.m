% test of wrapper for cost function
figure
p.Nmeas     = 8;
p.Nacro     = 32;
p.Nresidual = 1e3;
p.Nperm     = 1e2;
p.freq      = 3.9;
p.Amp       = 1.5;
p.noise     = 1.5;
p.Nbatch    = 10;

addpath('../utils_core')

clf
[t_unif,~]=getSamplingSchedules(p.Nmeas,0,0,0);
[I3,I4]=constructUtilMats(p);
N=5;
sc=1e-1;    
% wrap_getCostFun(t_unif,p,I3,I4)
tic
tiledlayout(2,1)
for kk=1:2
    if kk==1
        p.permMethod= 'slow';
    else
        p.permMethod= 'fast';
    end    
    nexttile(kk)
    disp(p.permMethod)
    % plot power only
    for ii=1:N
        [pwr,~,~]=getPowerBatch(t_unif,p,I3,I4);
        [~,ind]=max(mean(pwr,1));
        disp(ind)
        plot(mean(pwr,1),'-k')

        hold on
        ylim([0,1])
        drawnow
    end
end
toc
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