% example with bad variance but non-overlapping confidence intervals
tic
N=20000;
[unif_vals,nu_vals]=estimateBayesianFIMdet(5,5,1/3,'nreps',N);
fprintf('nu    mean:  %2.3f   std:  %2.3f  \n',mean(nu_vals),std(nu_vals)/sqrt(N)); % todo: add sample sdev
fprintf('unif  mean:  %2.3f   std:  %2.3f  \n\n',mean(unif_vals),std(unif_vals)/sqrt(N));
toc

function [mean_unif,mean_nu] = get_mean(NL,NR,tau)
N=20000;
[unif_vals,nu_vals]=estimateBayesianFIMdet(5,5,1/3,'nreps',N);
mean_unif=mean()
end