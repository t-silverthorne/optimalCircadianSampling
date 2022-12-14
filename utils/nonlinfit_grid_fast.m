function [best_fitresult,best_gof,ii_best,jj_best]  = nonlinfit_grid_fast(zts, Y, regularize)
% use linear regression to determine a good initial guess of the
% parameters, and then run nonlinear regression with this initial guess
dp=.1;
pergrid=.1:dp:1;

for ii=pergrid
    for jj=ii+dp:1
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
        if regularize
            if ii==min(pergrid) && jj==min(pergrid)+dp
                ii_best=ii;
                jj_best=jj;
                best_SSres=SSres;
            elseif SSres*(1-.5*(1-tanh(ii))*(1-tanh(jj))) ...
                < best_SSres*(1-.5*(1-tanh(ii_best))*(1-tanh(jj_best))) 
                ii_best=ii;
                jj_best=jj;
                best_SSres=SSres;
            end
        else
            if ii==min(pergrid) && jj==min(pergrid)+dp
                ii_best=ii;
                jj_best=jj;
                best_SSres=SSres;
            elseif SSres < best_SSres 
                ii_best=ii;
                jj_best=jj;
                best_SSres=SSres;
            end

        end
    end
end
[best_fitresult,best_gof]=nonlinfit(repmat(zts,size(Y,1),1), Y, ii_best,jj_best);

end

