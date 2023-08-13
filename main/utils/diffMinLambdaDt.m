function DlamDt = diffMinLambdaDt(mt,freq)
% derivative of smallest eigenvalue wrt measurement times 
X         = constructReducedX(mt,freq);
[W,D]     = eig(X'*X);
[~,emin]  = min(diag(D));
v         = W(:,emin);
v         = v/norm(v);

t3  = reshape(mt,[1,1,length(mt)]);
pf2 = 2*pi*freq;
dX  =pf2*[-2*cos(pf2*t3).*sin(pf2*t3)    cos(pf2*t3).^2 - sin(pf2*t3).^2;
      cos(pf2*t3).^2 - sin(pf2*t3).^2   2*cos(pf2*t3).*sin(pf2*t3)];

DlamDt = reshape(pagemtimes(pagemtimes(v',dX),v),[length(mt), 1]);
end
