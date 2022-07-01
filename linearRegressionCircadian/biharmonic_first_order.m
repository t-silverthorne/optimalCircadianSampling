%% First order iteration, no deletion $\alpha_s = (s+1)^{-1}
clear all
clf
mu=@(theta)  [1 cos(2*pi*theta) sin(2*pi*theta) cos(3.3*pi*theta) sin(3.3*pi*theta)]'*[1 cos(2*pi*theta) sin(2*pi*theta) cos(3.3*pi*theta) sin(3.3*pi*theta)];
psi=@(theta,theta_vec,w_vec) trace(mu(theta)*get_info_matrix(theta_vec,w_vec));


tiledlayout('flow')

theta_vec=[0];
w_vec=[1];
M_list=[]
for s=1:30
    xs=fminbnd(@(x) psi(x,theta_vec,w_vec),0,1);
    alphas=1/(s+1);
    w_vec=(1-alphas)*w_vec;
    w_vec(end+1)=alphas;
    theta_vec(end+1)=xs;
    M_list(end+1)=det(get_info_matrix(theta_vec,w_vec));
end
nexttile
polarscatter(theta_vec*2*pi,1,'.k')
nexttile
plot(1:numel(M_list),M_list)


theta_vec=[0];
w_vec=[1];
M_list=[];
psi=@(theta,theta_vec,w_vec) -5 + trace(mu(theta)*get_info_matrix(theta_vec,w_vec)/10);
for s=1:30
    xs_plus=fminbnd(@(x) psi(x,theta_vec,w_vec),0,1);
    psi_known=arrayfun( @(ind) psi(theta_vec(ind),theta_vec,w_vec),1:numel(theta_vec) );
    [psi_xs_minus,xs_minus_ind]=max(psi_known);
    xs_minus=theta_vec(xs_minus_ind);

    [psi(xs_plus,theta_vec,w_vec),-psi_xs_minus]
    [~,xs_ind]=min([psi(xs_plus,theta_vec,w_vec),-psi_xs_minus]);

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

nexttile
polarscatter(theta_vec*2*pi,1,'.k')
nexttile
plot(1:numel(M_list),M_list)

% %%
% disp('hi')
% theta_vec=[0];
% w_vec=[1];
% M_list=[];
% for s=1:10
%     disp(s)
%     xopt=simulannealbnd(@(x) wrap_psi_order_2(x,theta_vec,w_vec),[0.25 0.75 0.5],[0 0 0],[1 1 1]);
%     theta1=xopt(1);
%     theta2=xopt(2);
%     p1=xopt(3);
%     alphas=1/(s+1);
%     w_vec=(1-alphas)*w_vec;
%     w_vec(end+1)=alphas*p1;
%     w_vec(end+1)=alphas*(1-p1);
%     theta_vec(end+1)=theta1;
%     theta_vec(end+1)=theta2;
%     M_list(end+1)=det(get_info_matrix(theta_vec,w_vec));
% end
% nexttile
% polarscatter(theta_vec*2*pi,1,'.k')
% nexttile
% plot(1:numel(M_list),M_list)
% theta_vec
% w_vec


function M = get_info_matrix(theta_vec,w_vec)
% TODO: check normalization
syms theta
assume(theta,'real')
x=@(theta) [1 cos(2*pi*theta) sin(2*pi*theta) cos(3.3*pi*theta) sin(3.3*pi*theta)];

M=0;
for i=1:numel(theta_vec)
    M=M+w_vec(i)*x(theta_vec(i))'*x(theta_vec(i));
end
end

function psi_val = wrap_psi_order_2(x,theta_vec,w_vec)
mu=@(theta)  [1 cos(2*pi*theta) sin(2*pi*theta) cos(4*pi*theta) sin(4*pi*theta)]'*[1 cos(2*pi*theta) sin(2*pi*theta) cos(4*pi*theta) sin(4*pi*theta)];
psi_order_2=@(theta1,theta2,p1,theta_vec,w_vec) trace((p1*mu(theta1)+(1-p1)*mu(theta2))*get_info_matrix(theta_vec,w_vec));
psi_val=psi_order_2(x(1),x(2),x(3),theta_vec,w_vec);
end