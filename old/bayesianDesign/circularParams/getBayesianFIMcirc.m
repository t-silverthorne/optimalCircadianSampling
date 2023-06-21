function M = getBayesianFIMcirc(Ntimes,model)
% parameterize acrophase directly instead of linearizing model
syms x phi1 phi2 A1 A2 T1 T2
xv=sym('x',[1 Ntimes]);
assume([x phi1 phi2 A1 A2],'real')
switch model
    case 'cosinorOneFreq'
        f=A1.*cos(2*pi*x./T1-phi1);
        gradf=gradient(f,[A1 phi1 T1]);
    case 'cosinorTwoFreq'
        f=A1.*cos(2*pi*x./T1-phi1)+A2.*cos(2*pi*x./T2-phi2);
        gradf=gradient(f,[A1 phi1 T1 A2 phi2 T2]);
end
    
GF=subs(gradf,x,xv);
M=matlabFunction(GF*GF');
%GFfun=matlabFunction(GF);
end

