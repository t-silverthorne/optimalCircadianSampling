S=epsloc(:,:,1)
wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2))
%%
wvec(1).^2/sum(wvec.^2)*(sqrt(det(Sigm1))^(-1)*(2*pi)^(-d/2)*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,1)),pagemldivide(Sigm1,S-muvec(:,1))) )) + ...
                         wvec(2).^2/sum(wvec.^2)*(sqrt(det(Sigm2))^(-1)*(2*pi)^(-d/2)*...
                                           exp(-0.5*pagemtimes(pagetranspose(S-muvec(:,2)),pagemldivide(Sigm2,S-muvec(:,2))) ))