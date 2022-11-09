function [bestfit,lin_result] = multiStartRegression(zts, Y, model)
dp=.1;
pergrid=.1:dp:1;

switch model
    case 'cosinorOneFreq'
        for ii=pergrid
            x1=sin(2*pi*zts/ii);
            x2=cos(2*pi*zts/ii);
            x0=ones(1,numel(zts));
            X=[x0' x1' x2'];
            betas=(X'*X)\X'*Y'; 
            fits=(X*betas)';
            SSres=sum((fits-Y).^2,2);
            if ii==min(pergrid) || SSres < best_SSres
                ii_best=ii;
                best_SSres=SSres;
                lin_result.per1guess=ii_best;
                lin_result.a1guess=betas(2);
                lin_result.a2guess=betas(3);
            end

        end

end
bestfit=nonlinfit(zts, Y, model, lin_result);
end

