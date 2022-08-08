clear all
syms beta0 beta1 beta2 beta3 beta4
syms T1 T2
syms t
%%
gr=simplify(gradient(beta0 + beta1*sin(2*pi*t/T1)+beta2*cos(2*pi*t/T1)+beta3*sin(2*pi*t/T2)+beta4*cos(2*pi*t/T2),[beta0 beta1 beta2 beta3 beta4 T1 T2]));
matlabFunction(gr,'Vars',[t beta1 beta2 beta3 beta4 T1 T2])
fim_nonlinear=@(t,beta1,beta2,beta3,beta4,T1,T2)reshape([1.0,sin((t.*pi.*2.0)./T1),cos((t.*pi.*2.0)./T1),sin((t.*pi.*2.0)./T2),cos((t.*pi.*2.0)./T2),1.0./T1.^2.*t.*pi.*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0,1.0./T2.^2.*t.*pi.*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0,sin((pi.*conj(t).*2.0)./conj(T1)),sin((pi.*conj(t).*2.0)./conj(T1)).*sin((t.*pi.*2.0)./T1),sin((pi.*conj(t).*2.0)./conj(T1)).*cos((t.*pi.*2.0)./T1),sin((pi.*conj(t).*2.0)./conj(T1)).*sin((t.*pi.*2.0)./T2),sin((pi.*conj(t).*2.0)./conj(T1)).*cos((t.*pi.*2.0)./T2),1.0./T1.^2.*t.*pi.*sin((pi.*conj(t).*2.0)./conj(T1)).*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0,1.0./T2.^2.*t.*pi.*sin((pi.*conj(t).*2.0)./conj(T1)).*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0,cos((pi.*conj(t).*2.0)./conj(T1)),cos((pi.*conj(t).*2.0)./conj(T1)).*sin((t.*pi.*2.0)./T1),cos((pi.*conj(t).*2.0)./conj(T1)).*cos((t.*pi.*2.0)./T1),cos((pi.*conj(t).*2.0)./conj(T1)).*sin((t.*pi.*2.0)./T2),cos((pi.*conj(t).*2.0)./conj(T1)).*cos((t.*pi.*2.0)./T2),1.0./T1.^2.*t.*pi.*cos((pi.*conj(t).*2.0)./conj(T1)).*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0,1.0./T2.^2.*t.*pi.*cos((pi.*conj(t).*2.0)./conj(T1)).*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0,sin((pi.*conj(t).*2.0)./conj(T2)),sin((pi.*conj(t).*2.0)./conj(T2)).*sin((t.*pi.*2.0)./T1),sin((pi.*conj(t).*2.0)./conj(T2)).*cos((t.*pi.*2.0)./T1),sin((pi.*conj(t).*2.0)./conj(T2)).*sin((t.*pi.*2.0)./T2),sin((pi.*conj(t).*2.0)./conj(T2)).*cos((t.*pi.*2.0)./T2),1.0./T1.^2.*t.*pi.*sin((pi.*conj(t).*2.0)./conj(T2)).*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0,1.0./T2.^2.*t.*pi.*sin((pi.*conj(t).*2.0)./conj(T2)).*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0,cos((pi.*conj(t).*2.0)./conj(T2)),cos((pi.*conj(t).*2.0)./conj(T2)).*sin((t.*pi.*2.0)./T1),cos((pi.*conj(t).*2.0)./conj(T2)).*cos((t.*pi.*2.0)./T1),cos((pi.*conj(t).*2.0)./conj(T2)).*sin((t.*pi.*2.0)./T2),cos((pi.*conj(t).*2.0)./conj(T2)).*cos((t.*pi.*2.0)./T2),1.0./T1.^2.*t.*pi.*cos((pi.*conj(t).*2.0)./conj(T2)).*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0,1.0./T2.^2.*t.*pi.*cos((pi.*conj(t).*2.0)./conj(T2)).*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0,pi.*1.0./conj(T1).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-2.0,pi.*sin((t.*pi.*2.0)./T1).*1.0./conj(T1).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-2.0,pi.*cos((t.*pi.*2.0)./T1).*1.0./conj(T1).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-2.0,pi.*sin((t.*pi.*2.0)./T2).*1.0./conj(T1).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-2.0,pi.*cos((t.*pi.*2.0)./T2).*1.0./conj(T1).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-2.0,1.0./T1.^2.*t.*pi.^2.*1.0./conj(T1).^2.*conj(t).*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-4.0,1.0./T2.^2.*t.*pi.^2.*1.0./conj(T1).^2.*conj(t).*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*(cos((pi.*conj(t).*2.0)./conj(T1)).*conj(beta1)-sin((pi.*conj(t).*2.0)./conj(T1)).*conj(beta2)).*-4.0,pi.*1.0./conj(T2).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-2.0,pi.*sin((t.*pi.*2.0)./T1).*1.0./conj(T2).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-2.0,pi.*cos((t.*pi.*2.0)./T1).*1.0./conj(T2).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-2.0,pi.*sin((t.*pi.*2.0)./T2).*1.0./conj(T2).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-2.0,pi.*cos((t.*pi.*2.0)./T2).*1.0./conj(T2).^2.*conj(t).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-2.0,1.0./T1.^2.*t.*pi.^2.*1.0./conj(T2).^2.*conj(t).*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-4.0,1.0./T2.^2.*t.*pi.^2.*1.0./conj(T2).^2.*conj(t).*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*(cos((pi.*conj(t).*2.0)./conj(T2)).*conj(beta3)-sin((pi.*conj(t).*2.0)./conj(T2)).*conj(beta4)).*-4.0],[7,7]);