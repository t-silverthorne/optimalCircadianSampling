clear
%% Illustrates non-monotonic convergence
k=6;
n0=1;
theta_vec=rand(n0,1);
w_vec=ones(n0,1)/numel(theta_vec);

for i=1:50
    Mnow=get_info_matrix_k(theta_vec,w_vec,k);
    f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];
    xnew=fminbnd(@(x) -trace(Mnow\(f(x)'*f(x))),0,1);
    theta_vec=[theta_vec; xnew];
    w_vec=ones(n0+i,1)/numel(theta_vec);
    semilogy(i,log(det(get_info_matrix_k(theta_vec,w_vec,k))),'ok')
    hold on
end
xline(k*(k+1)/2)
hold off

function M = get_info_matrix_k(theta_vec,w_vec,k)
syms theta
assume(theta,'real')
f=@(theta) [1 cos(2*pi*theta*(1:k)) sin(2*pi*theta*(1:k))];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*f(theta_vec(i))'*f(theta_vec(i));
end
end
