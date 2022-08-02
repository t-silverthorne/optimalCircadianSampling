clear
X=open('Nsamp_1_nreps_100_per1_12_per2_4.mat')
%%
sse_unif=NaN(numel(X.gof_unif),1);
sse_nu=NaN(numel(X.gof_nu),1);

for i=1:numel(sse_unif)
    sse_unif(i)=X.gof_unif{i}.sse;
    sse_nu(i)=X.gof_nu{i}.sse;
end
close all
histogram(sse_unif)
hold on
histogram(sse_nu)
legend({'unif','nu'})
%%
per_unif=NaN(numel(X.gof_unif),1);
per_nu=NaN(numel(X.gof_nu),1);

for i=1:numel(per_unif)
    per_unif(i)=max(X.res_unif{i}.per1,X.res_unif{i}.per2);
    per_nu(i)=max(X.res_nu{i}.per1,X.res_nu{i}.per2);
end
close all
histogram(per_unif)
hold on
histogram(per_nu)
legend({'unif','nu'})