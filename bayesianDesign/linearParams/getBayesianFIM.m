function M = getBayesianFIM(Ntimes)
% Fisher information matrix for a biharmonic cosinor model
% calculated symbolically and returned as a MATLAB function
syms x beta1 beta2 beta3 beta4 beta5 beta6
xv=sym('x',[1 Ntimes]);
assume([x beta1 beta2 beta3 beta4 beta5 beta6],'real')
f=beta1.*sin(2*pi*x./beta5)+beta2.*cos(2*pi*x./beta5) + ...
                     beta3.*sin(2*pi*x./beta6)+beta4.*cos(2*pi*x./beta6);
gradf=gradient(f,[beta1 beta2 beta3 beta4 beta5 beta6]);
GF=subs(gradf,x,xv);
GFfun=matlabFunction(GF);
M=matlabFunction(GF*GF');
end

