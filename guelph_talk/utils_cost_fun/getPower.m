function [pwr_est,amp_est,phi_est,bin_mat] = getPower(t,p)
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

YI        = getPermDataAllMethods(p,Y);

betas     = pagemldivide(X'*X,pagemtimes(X',pagetranspose(YI)));
fits      = pagetranspose(pagemtimes(X,betas));
SSres     = sqrt(sum((fits-YI).^2,2));
bin_mat   = SSres<SSres_obs;
pwr_est   = sum(sum(bin_mat,3)/p.Nperm<.05)/p.Nresidual;
end
