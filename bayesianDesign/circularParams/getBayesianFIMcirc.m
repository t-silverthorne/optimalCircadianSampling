function M = getBayesianFIMcirc(Ntimes,nharm)
% parameterize acrophase directly instead of linearizing model
syms x phi1 phi2 A1 A2 T1 T2
xv=sym('x',[1 Ntimes]);
assume([x phi1 phi2 A1 A2],'real')
if nharm==1
    f=A1.*cos(2*pi*x./T1-phi1);
    gradf=gradient(f,[A1 phi1 T1]);
elseif nharm==2
    f=A1.*cos(2*pi*x./T1-phi1)+A2.*cos(2*pi*x./T2-phi2);
    gradf=gradient(f,[A1 phi1 T1 A2 phi2 T2]);
else
    f=NaN;gradf=NaN;disp('Error: can only handle one or two frequencies')
end
    
GF=subs(gradf,x,xv);
GFfun=matlabFunction(GF);
M=matlabFunction(GF*GF');
end

