function M = get_FIM_biharmonic_nonlin(theta_vec,w_vec,param)
f_nonlinear=@(t,beta1,beta2,beta3,beta4,T1,T2)[1.0;sin((t.*pi.*2.0)./T1);cos((t.*pi.*2.0)./T1);sin((t.*pi.*2.0)./T2);cos((t.*pi.*2.0)./T2);1.0./T1.^2.*t.*pi.*(beta2.*sin((t.*pi.*2.0)./T1)-beta1.*cos((t.*pi.*2.0)./T1)).*2.0;1.0./T2.^2.*t.*pi.*(beta4.*sin((t.*pi.*2.0)./T2)-beta3.*cos((t.*pi.*2.0)./T2)).*2.0];
x=@(theta) f_nonlinear(theta,param.beta1,param.beta2, ...
           param.beta3,param.beta4,param.T1,param.T2);
M=0;
for ii=1:numel(theta_vec)
    M=M+w_vec(ii)*x(theta_vec(ii))*x(theta_vec(ii))';
end
end