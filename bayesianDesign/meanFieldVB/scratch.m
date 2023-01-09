
% prior is product of Gaussians
prior=@(A1,B1,T1) exp(-0.5*([A1;B1;T1] - muvec_prior)'*(Sigmat_prior\([A1;B1;T1] - muvec_prior)))*(2*pi)^(-d/2)/sqrt(det(Sigmat_prior));
lhood=@(tobs,yobs,A1,B1,T1) exp(-0.5*sum(( A1*cos(pi*fmax*(1+tanh(T1)*tobs))+B1*sin(pi*fmax*(1+tanh(T1))*tobs) -yobs).^2))/sqrt(2*pi);
% vectorized version of prior (capable of using higher dimensional arrays)
priorvec=@(S) exp(pagemtimes(-0.5*pagetranspose(S - muvec_prior), ...
                             (pagemldivide(Sigmat_prior,(S - muvec_prior)))))*...
                  (2*pi)^(-d/2)/sqrt(det(Sigmat_prior));

eta=@(S,tobs) S(1,:,:).*cos(pi*fmax*(1+tanh(S(3,:,:)).*tobs))+...
    S(2,:,:).*sin(pi*fmax*(1+tanh(S(3,:,:))).*tobs);
lhoodvec=@(S,tobs,yobs) exp(-0.5*sum( (eta(S,tobs)-yobs).^2,1))/sqrt(2*pi);
