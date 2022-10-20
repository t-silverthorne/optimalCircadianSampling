function M = getMultiNormFIM(N,method)
mu=sym('mu',[N 1]);
sig=sym('s',[N N]);
xv=sym('xv',[N 1]);
%%
f=-(xv-mu)'*inv(sig)*(xv-mu)/sqrt((2*pi)^N*det(sig))
grad(f,mu)
end

