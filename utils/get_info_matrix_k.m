function M = get_info_matrix_k(theta_vec,w_vec,k)
f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*f(theta_vec(i))'*f(theta_vec(i));
end
end
