function M = get_info_matrix(theta_vec,w_vec)
% Information matrix for one-frequency harmonic regression
syms theta
assume(theta,'real')
x=@(theta) [1 cos(2*pi*theta) sin(2*pi*theta)];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*x(theta_vec(i))'*x(theta_vec(i));
end
end
