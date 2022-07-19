% ignore problem structure, just do naive approach to minimize det(M)
clear
mu=@(theta) [1 cos(2*pi*theta) sin(2*pi*theta)]'*[1 cos(2*pi*theta) sin(2*pi*theta)];
psi=@(theta,theta_vec,w_vec) trace(mu(theta)*get_info_matrix(theta_vec,w_vec));

num_points=9;
theta_vec=2*pi*rand(num_points,1);
w_vec=1/num_points*ones(num_points,1);%diff( [0 sort(rand(1,num_points-1)) 1])'

det(get_info_matrix(theta_vec,w_vec))
thetas_best=linspace(0,2*pi,num_points);
thetas_best=thetas_best(1:end-1);
det(get_info_matrix(thetas_best, ones(numel(thetas_best),1)*1/numel(thetas_best)))
%%

det(get_info_matrix([0 .25 .5 .75]',ones(4,1)/4))
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

function M = get_info_matrix2(theta_vec,w_vec)
% TODO: check normalization
syms theta
assume(theta,'real')
x=@(theta) [ones(numel(theta),1) cos(2*pi*theta) sin(2*pi*theta)];
VSelect=@(V,IDX) V(:,IDX);
M=NaN(3,3);
for i=1:numel(x(0))
    for j=1:numel(x(0))
        M(i,j)=VSelect(x(theta_vec),i)'*( VSelect(x(theta_vec),j).*w_vec );
    end
end
end

