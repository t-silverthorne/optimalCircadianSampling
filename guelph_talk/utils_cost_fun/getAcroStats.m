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
acro_bias       = mean(abs(exp(1j*acro_mat)- exp(1j*phivec)));
acro_var        = abs(mean(exp(1j*acro_mat),1));
end

