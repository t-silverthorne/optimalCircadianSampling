function [amp_bias,amp_var] = getAmpStats(amp_mat,A0)
%%%%%%%%%%%%%%
%GETAMPSTATS compute statistics of acrophase estimates
% INPUT: 
%   amp_mat   matrix of amplitude estimates
%   A0        true amplitude of cosinor model
% OUTPUT:
%   amp_bias  E[\hat{A} - A0]   
%   amp_var   E[(\hat{A}-A0)^2]
%%%%%%%%%%%%%%
amp_mean=mean(amp_mat,1);
amp_bias=abs(amp_mean-A0);
amp_var =var(amp_mat);
end

