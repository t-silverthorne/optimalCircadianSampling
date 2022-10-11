% example with bad variance but non-overlapping confidence intervals
tic
[unif_vals,nu_vals]=estimateBayesianFIMdet(5,5,1/3,'nreps',30000);
[~,~,nuci]=normfit(nu_vals);
[~,~,unifci]=normfit(unif_vals);
fprintf('nu    mean:  %2.3f   std:  %2.3f  CI: [%2.3f,%2.3f]\n',mean(nu_vals),std(nu_vals),nuci(1),nuci(2)); % todo: add sample sdev
fprintf('unif  mean:  %2.3f   std:  %2.3f  CI: [%2.3f,%2.3f]\n\n',mean(unif_vals),std(unif_vals),unifci(1),unifci(2));
toc