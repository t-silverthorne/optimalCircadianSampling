function [best_fitresult,best_gof] = nonlinfit_grid(zts, Xdat)
pergrid=1:1:24;

best_fitresult=NaN;
best_gof=NaN;
for ii=pergrid
    for jj=ii+1:24
        [fitresult,gof]=nonlinfit(zts, Xdat, ii, jj);
        if ii==1 && jj==2
            best_fitresult=fitresult;
            best_gof=gof;
        elseif gof.adjrsquare > best_gof.adjrsquare
            best_fitresult=fitresult;
            best_gof=gof;
        end
    end
end

end

