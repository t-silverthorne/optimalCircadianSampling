function [best_fitresult,best_gof,ii_best,jj_best]  = nonlinfit_grid_fast(zts, Y)
pergrid=1:1:24;

for ii=pergrid
    for jj=ii+1:24
        x1=sin(2*pi*zts/ii);
        x2=cos(2*pi*zts/ii);
        x0=ones(1,numel(zts));
        x3=sin(2*pi*zts/jj);
        x4=cos(2*pi*zts/jj);
        X= [x0' x1' x2' x3' x4'];
        % do linear regression 
        betas=(X'*X)\X'*Y'; 
        fits=(X*betas)';
        SSres=sum((fits-Y).^2,2);
        if ii==1 && jj==2
            ii_best=ii;
            jj_best=jj;
            best_SSres=SSres;
        elseif SSres*(1-.5*(1-tanh(ii))*(1-tanh(jj))) < best_SSres
            ii_best=ii;
            jj_best=jj;
            best_SSres=SSres;
        end
    end
end
[best_fitresult,best_gof]=nonlinfit(zts, Y, ii_best,jj_best);

end

