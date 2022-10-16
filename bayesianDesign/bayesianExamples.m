% example with bad variance but non-overlapping confidence intervals
close all
%% 1D plot
tiledlayout(1,3)
for jj=[1/4 1/3 1/2]
    nexttile
    for i=3:10
        [unif_vals,sd_unif,nu_vals,sd_nu]=estimateBayesianFIMdet(i,2+i,jj,'sdev');
        mean_unif=mean(unif_vals)
        mean_nu=mean(nu_vals)
        semilogy(i,mean_unif,'.k','MarkerSize',20)
        hold on
        semilogy(i,mean_unif+sd_unif,'*k','MarkerSize',2)
        semilogy(i,mean_unif-sd_unif,'*k','MarkerSize',2)
    
        semilogy(i,mean_nu,'.b','MarkerSize',20)
        semilogy(i,mean_nu+sd_nu,'*b','MarkerSize',2)
        semilogy(i,mean_nu-sd_nu,'*b','MarkerSize',2)
        drawnow
    end
    hold off
    ylim([20 40])
end

%% heatmap
nii=8;
njj=8;
unif=NaN(nii,njj);
nu=NaN(nii,njj);
for ii=1:nii
    disp(ii)
    for jj=1:njj
        [unif_vals,~,nu_vals,~]=estimateBayesianFIMdet(i,2+i,jj,'sdev');
        unif(ii,jj)=mean(unif_vals);
        nu(ii,jj)=mean(nu_vals);
    end
end
%%
close all
tiledlayout(1,2)
nexttile
heatmap(1:njj,1:nii,real(unif))
nexttile
heatmap(1:njj,1:nii,real(nu))


%%
%fprintf('nu    mean:  %2.3f   std:  %2.3f  \n',mean(nu_vals),std(nu_vals)/sqrt(N)); % todo: add sample sdev
%fprintf('unif  mean:  %2.3f   std:  %2.3f  \n\n',mean(unif_vals),std(unif_vals)/sqrt(N));
%toc

function [mean_unif,mean_nu] = get_mean(NL,NR,tau)
N=20000;
[unif_vals,nu_vals]=estimateBayesianFIMdet(5,5,1/3,'nreps',N);
mean_unif=mean()
end