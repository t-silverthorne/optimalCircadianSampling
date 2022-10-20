function M = getCosinorTwoFreqFIM(Ntimes,method)
% options for method: hessian, grad
if nargin<2
    method='hessian';
end
% parameterize acrophase directly instead of linearizing model
syms x phi1 A1 T1 phi2 A2 T2
xv=sym('x',[1 Ntimes]);
assume([x phi1 A1],'real')
f=A1.*cos(2*pi*x./T1-phi1);
switch method
    case 'hessian'
        f=A1.*cos(2*pi*x./T1-phi1)+A2.*cos(2*pi*x./T2-phi2);
        hessf=hessian(f,[A1 phi1 T1 A2 phi2 T2]);
        M=0;
        for ii=1:numel(xv)
            M=M+subs(hessf,x,xv(ii));
        end
        M=matlabFunction(M);
    case 'grad'
        gradf=gradient(f,[A1 phi1 T1 A2 phi2 T2]);
        GF=subs(gradf,x,xv);
        M=matlabFunction(GF*GF');
    case 'gradsum'
        gradf=gradient(f,[A1 phi1 T1 A2 phi2 T2]);
        M=0;
        for ii=1:numel(xv)
            M=M+subs(gradf,x,xv(ii))*(subs(gradf,x,xv(ii))');
        end
        M=matlabFunction(M);
end

end

