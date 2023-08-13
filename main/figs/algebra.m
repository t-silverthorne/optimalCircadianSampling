N=2;
A=rand(N,N);
A=A'*A;
B= rand(N,N);
C=rand(N,N);
B=B'*B;C=C'*C;
eA = eig(A);
c1 = max(eA);
c2 = randsample(1:(N-1),1);
c2 = eA(c2);


chol( B*pinv(c1*eye(N)-A)*C + C'*pinv(c1*eye(N)-A)*B' - ( B*pinv(c2*eye(N)-A)*C + C'*pinv(c2*eye(N)-A)*B' ) )
%%
norm(pinv(c1*eye(N)-A)) < norm(pinv(c2*eye(N)-A))

%%
eA
norm(pinv(A))
1/min(eA)
%%
N=2
X=constructReducedX(rand(1,N),rand)
trace(inv(X'*X))
N/det(X'*X)
N*det(inv(X'*X))

eig(inv(X'*X))
1/N

%% Verify acrophase estimate
N=randsample(1:1e2,1);
u=10^rand*rand(N,1);
v=10^rand*rand(N,1);

e=ones(N,1)/sqrt(N);
1-v'*u
(e-v)'*(e+u)- e'*u+e'*v