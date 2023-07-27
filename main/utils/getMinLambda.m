function lmin=getMinLambda(freq,mt)
% worst case lambda for all acrophases is min(eig(X'*X))
X=constructReducedX(mt,freq);
lmin = min(eig(X'*X));
end
