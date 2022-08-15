% old code, use the implementation in the linearRegressionCircadian/clean/ directory

%% First order iteration, no deletion $\alpha_s = (s+1)^{-1}
mu=@(theta) [1 cos(2*pi*theta) sin(2*pi*theta)]'*[1 cos(2*pi*theta) sin(2*pi*theta)];
psi=@(theta,theta_vec,w_vec) trace(mu(theta)*get_info_matrix(theta_vec,w_vec));

theta_vec=[0];
w_vec=[1];
M_list=[]
for s=1:3
    xs=fminbnd(@(x) psi(x,theta_vec,w_vec),0,1);
    alphas=1/(s+1);
    w_vec=(1-alphas)*w_vec;
    w_vec(end+1)=alphas;
    theta_vec(end+1)=xs;
    M_list(end+1)=det(get_info_matrix(theta_vec,w_vec));
end
theta_vec
w_vec
tiledlayout('flow')
nexttile
polarscatter(theta_vec*2*pi,1,'.k')
nexttile
plot(1:numel(M_list),M_list)

theta_vec=linspace(0,1,5);
theta_vec=theta_vec(1:end-1);
w_vec=ones(1,numel(theta_vec))/numel(theta_vec);
d_opt=det(get_info_matrix(theta_vec,w_vec))
yline(d_opt)
%% First order iteration, now with deletion
%theta_vec=linspace(0,1,6)+.05*randn(1,6);
%theta_vec=theta_vec(1:end-1);
theta_vec=[0];
w_vec=[1];
%w_vec=ones(1,numel(theta_vec))/numel(theta_vec);
M_list=[]
for s=1:20
    xs_plus=fminbnd(@(x) psi(x,theta_vec,w_vec),0,1);
    psi_known=arrayfun( @(ind) psi(theta_vec(ind),theta_vec,w_vec),1:numel(theta_vec) );
    [psi_xs_minus,xs_minus_ind]=max(psi_known);
    xs_minus=theta_vec(xs_minus_ind);
    if s<5
        xs_ind=1;
    else
        [~,xs_ind]=min([psi(xs_plus,theta_vec,w_vec),-psi_xs_minus]);
    end
    if xs_ind==1
        xs=xs_plus;
        alphas=1/(s+1);
        w_vec=(1-alphas)*w_vec;
        w_vec(end+1)=alphas;
        theta_vec(end+1)=xs;
    else
        xs=xs_minus;
        alphas=-min(1/(s+1), w_vec(xs_minus_ind)/(1-w_vec(xs_minus_ind)) );
        w_vec=(1-alphas)*w_vec;
        w_vec(xs_minus_ind)=w_vec(xs_minus_ind)+alphas;
    end
    M_list(end+1)=det(get_info_matrix(theta_vec,w_vec));
end
theta_vec
w_vec
tiledlayout('flow')
nexttile
polarscatter(theta_vec*2*pi,1,'.k')
nexttile
plot(1:numel(M_list),M_list)


function M = get_info_matrix(theta_vec,w_vec)
% TODO: check normalization
syms theta
assume(theta,'real')
x=@(theta) [1 cos(2*pi*theta) sin(2*pi*theta)];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*x(theta_vec(i))'*x(theta_vec(i));
end
end