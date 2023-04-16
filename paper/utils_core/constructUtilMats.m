function [X,I3,I4] = constructUtilMats(t,p)
%%%%%%%%%%%%%%
%CONSTRUCTUTILMATS Construct matrices used throughout the program, X is
% used in linear regression, and I3,I4 are used in getPermutedData().
% INPUT:   struct p is the parameters
%   t               measurement times
%   p.Nmeas         number of measurements in each experiment
%   p.Nacro         number of acrophases to simulate
%   p.Nresidual     number of simulated experiment (for each acrophase)
%   p.Nperm         number of permutations for hypothesis testing
%   p.freq          frequency of cosinor model
% OUTPUT:
%   X               design matrix used in harmonic regression
%   I3              index matrix used in sub2ind() permutation operations
%   I4              index matrix used in sub2ind() permutation operations
%%%%%%%%%%%%%%

% build design matrix X
x0=ones(1,length(t));
x1=sin(2*pi*p.freq*t);
x2=cos(2*pi*p.freq*t);
X= [x0' x1' x2'];

% build I3 and I4, bookeeping matrices for permutations
N1=p.Nresidual; N2=p.Nmeas; N3=p.Nperm; N4=p.Nacro;
I3=NaN(N1,N2,N3);
I4=NaN(N1,N2,N3,N4);
for ii=1:N3
    I3(:,:,ii)=ii*ones(N1,N2);
end
for ii=1:N4
    I4(:,:,:,ii)=ii*ones(N1,N2,N3);
end
end

