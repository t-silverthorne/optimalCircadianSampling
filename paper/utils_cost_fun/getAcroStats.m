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
acro_mean  = angle(mean(exp(1j*acro_mat),1));
acro_bias  = min([mod(acro_mean-phivec,2*pi);mod(phivec-acro_mean,2*pi)],[],1);
acro_var   = arrayfun(@norm,acro_mean.*exp(-1j*phivec));
end

