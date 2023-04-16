function YI = getPermutedData(Y,R,I3,I4)
%%%%%%%%%%%%%%
%GETPERMUTEDDATA act on simulated data Y with permutations stored in R
% INPUT:   
%   Y          simulated data (size Nresidual x Nmeas x 1 x Nacro)
%   R          permutations   (size Nresidual x Nmeas x Nperm x Nacro)
%   I3         bookkeeping matrix for sub2ind, generated in constructUtilMats
%   I4         bookkeeping matrix for sub2ind, generated in constructUtilMats
% OUTPUT:
%   YI         permuted data (size Nresidual x Nmeas x Nperm x Nacro)
%%%%%%%%%%%%%%
N1=size(R,1);N2=size(R,2);N3=size(R,3);N4=size(R,4);
Y=repmat(Y,1,1,N3,1);
YI=Y(sub2ind(size(Y),repmat((1:N1)',1,N2,N3,N4),R,repmat(I3,1,1,1,N4),I4));
end
