function [unif_vals,nu_vals] = getSimulatedData(M,nreps,mt_unif,mt_nu,run_unif,run_nu)
% generate simulated data from biharmonic cosinor model and then return a
% vector of estimates for log(det(M(xi,beta)))
% 
% the calculation can be run with either a pre-determined number of samples
% or until some variance threshold has been passed 

% inputs
unif_vals=NaN; nu_vals=NaN;
if nargin < 5
    run_unif=true; run_nu=true;
end


% generate data
beta=rand(nreps,4);% first four parameters are unif [10,11]
beta(:,1:4)=normrnd(10,1,nreps,4);
accepted=NaN(nreps,2); 
% final two need rejection sampling
ind=1;
while ind<nreps
    samp=rand(floor(nreps*1.1),2);
    acceptedloc=samp(abs(samp(:,1)-samp(:,2))>.2,:);
    accepted(ind:ind+size(acceptedloc,1)-1,:)=acceptedloc;
    ind=ind+size(acceptedloc,1);
end
beta(:,5:6)=accepted(1:nreps,:);
% beta=NaN(nreps,6);% first four parameters are unif [0,1]
% beta(:,1:4)=normrnd(10,1,nreps,4);

% get log(det(M)) for each parameter set
beta=num2cell(beta);
mtc_unif=num2cell(mt_unif);
mtc_nu=num2cell(mt_nu);
if run_unif 
    unif_vals=arrayfun(@(ind) log(det(M(beta{ind,:},mtc_unif{:}))),1:nreps);
end
if run_nu
    nu_vals=arrayfun(@(ind) log(det(M(beta{ind,:},mtc_nu{:}))),1:nreps);
end
end

