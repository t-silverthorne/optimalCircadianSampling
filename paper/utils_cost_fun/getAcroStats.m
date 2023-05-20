function [acro_bias,acro_var] = getAcroStats(acro_mat,phivec)
%%%%%%%%%%%%%%
%GETACROSTATS compute statistics of acrophase estimates
% INPUT: 
%   acro_mat   matrix of acrophase estimates
%   phi0       true acrophase of cosinor model
% OUTPUT:
%   acro_bias  expected value of arclength distance between estimator and phi0
%   acro_var   complex modulus of E[ e^{i (phi_k - phi_0)}]
%%%%%%%%%%%%%%
zacro_mean  = mean(exp(1j*acro_mat),1);
acro_mean_angle =angle(zacro_mean);
acro_bias  = abs(exp(1j*acro_mean_angle)- exp(1j*phi_vec));
acro_var   = arrayfun(@norm,zacro_mean.*exp(-1j*phivec));
end

