function [unif_vals,nu_vals] = getSimulatedDataCirc(M,nreps,mt_unif,mt_nu,run_unif,run_nu)
%GETSIMULATEDDATACIRC Summary of this function goes here
%   Detailed explanation goes here
% inputs
unif_vals=NaN; nu_vals=NaN;
if nargin < 5
    run_unif=true; run_nu=true;
end

% generate data [A1 phi1 T1 A2 phi2 T2]
beta=NaN(nreps,6); % parameter matrix
beta(:,[1 4])=10*rand(nreps,2); % get amplitudes
beta(:,[2 5])= 2*pi*rand(nreps,2); % get acrophases
accepted=NaN(nreps,2); 
% final two need rejection sampling
ind=1;
while ind<nreps
    samp=rand(floor(nreps*1.1),2);
    acceptedloc=samp(abs(samp(:,1)-samp(:,2))>.2,:);
    accepted(ind:ind+size(acceptedloc,1)-1,:)=acceptedloc;
    ind=ind+size(acceptedloc,1);
end
beta(:,[3 6])=accepted(1:nreps,:); % get periods

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

