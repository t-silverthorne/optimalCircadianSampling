function M = getBayesianFIMcirc(Ntimes,model)
% parameterize acrophase directly instead of linearizing model
syms x phi1 phi2 A1 A2 f1 f2
xv=sym('x',[1 Ntimes]);
assume([x phi1 phi2 A1 A2 f1 f2],'real')
switch model
    case 'cosinorOneFreq'
        f=A1.*cos(2*pi*x.*f1-phi1);
        gradf=gradient(f,[A1 phi1 f1]);
        
        GF=subs(gradf,x,xv);
        M=matlabFunction(GF*GF','Vars',[A1 phi1 f1 xv]);
    case 'cosinorTwoFreq'
        f=A1.*cos(2*pi*x.*f1-phi1)+A2.*cos(2*pi*x.*f2-phi2);
        gradf=gradient(f,[A1 phi1 f1 A2 phi2 f2]);

        GF=subs(gradf,x,xv);
        M=matlabFunction(GF*GF','Vars',[A1 phi1 f1 A2 phi2 f2 xv]);
end

%GFfun=matlabFunction(GF);
end

