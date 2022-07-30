function M = get_memory_info_matrix(theta_vec,w_vec,kernel_opt)
% TODO: check normalization
syms theta
assume(theta,'real')
x=@(theta) [1 cos(2*pi*theta) sin(2*pi*theta)];
sigma0=1;
gamma0=.5;

switch kernel_opt
    case 'exp'
        alpha = @(x,xp) exp(-abs(.25*(x-xp)));
    case 'tanhsq'
        alpha = @(x,xp) 1-tanh(abs(.25*(x-xp))).^2;
end

M=0;
for i=1:numel(theta_vec)
    ss=sigma0*(1+ gamma0*sum(w_vec(1:i-1).*alpha(theta_vec(i),theta_vec(1:i-1))));
    M=M+w_vec(i)*x(theta_vec(i))'*x(theta_vec(i))/ss^2;
end
end
