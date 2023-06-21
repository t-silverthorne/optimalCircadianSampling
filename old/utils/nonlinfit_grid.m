function [best_fitresult,best_gof,ii_best,jj_best] = nonlinfit_grid(zts, Xdat)
% nonlinear regression with multiple initial guesses of the periods
% NOTE: this is orders of magnitude slower than nonlinfit_grid_fast.m
pergrid=1:1:24;

best_fitresult=NaN;
best_gof=NaN;
for ii=pergrid
    for jj=ii+1:24
        [fitresult,gof]=nonlinfit(zts, Xdat, ii, jj);
        if ii==1 && jj==2
            ii_best=ii;
            jj_best=jj;
            best_fitresult=fitresult;
            best_gof=gof;
        elseif gof.adjrsquare > best_gof.adjrsquare
            ii_best=ii;
            jj_best=jj;
            best_fitresult=fitresult;
            best_gof=gof;
        end
    end
end

end

