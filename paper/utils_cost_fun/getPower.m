function [pwr_est,amp_est,phi_est] = getPower(t,p,I3,I4)
%%%%%%%%%%%%%%
%GETPOWER compute power as a function of acrophase
% INPUT:   
%   t          vector of measurement times
%   p          param struct, storing size of simulation and osc features
%   I3         bookkeeping matrix for sub2ind, generated in constructUtilMats
%   I4         bookkeeping matrix for sub2ind, generated in constructUtilMats
% OUTPUT:
%   pwr_est    estimate of power as a function of acrophase
%   amp_est    amplitude estimates for each acrophase
%   phi_est    acrophase estimates for each acrophase
%%%%%%%%%%%%%%
Y=getSimulatedData(t,p);
X=constructX(t,p);

betas_obs = pagemldivide(X'*X,pagemtimes(X',pagetranspose(Y))); % observed error
phi_est   = atan2(betas_obs(2,:,:,:),betas_obs(3,:,:,:));
amp_est   = sqrt(betas_obs(2,:,:,:).^2 + betas_obs(3,:,:,:).^2);
fits_obs  = pagetranspose(pagemtimes(X,betas_obs));
SSres_obs = sqrt(sum((fits_obs-Y).^2,2));
clear betas_obs fits_obs

R  = getPermutations(p.Nresidual,p.Nmeas,p.Nperm,p.Nacro);
YI = getPermutedData(Y,R,I3,I4);

betas   = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits    = pagetranspose(pagemtimes(X,betas));
SSres   = sqrt(sum((fits-YI).^2,2));
pwr_est = sum(sum(SSres>SSres_obs,3)/p.Nperm>.95)/p.Nresidual;
end
