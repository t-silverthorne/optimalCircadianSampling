function [period_error]=get_mean_period_error(zts,Y,per1,per2)
% for fitmethod2
modelfun = @(beta,x) beta(1).*sin(2*pi*x./beta(5))+beta(2).*cos(2*pi*x./beta(5)) + ...
                     beta(3).*sin(2*pi*x./beta(6))+beta(4).*cos(2*pi*x./beta(6)) ;
fitmethod = 2;
dp=.1;
pergrid=.1:dp:1;
for ii=pergrid
    for jj=ii+dp:dp:1
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
        SSresavg=mean(SSres);
        
        if ii==min(pergrid) && jj==min(pergrid)+dp
            ii_best=ii;
            jj_best=jj;
            best_SSresavg=SSresavg;
        elseif SSresavg < best_SSresavg 
            ii_best=ii;
            jj_best=jj;
            best_SSresavg=SSresavg;
        end
    end
end
perSmall_vals=NaN(1,size(Y,1));
perBig_vals=NaN(1,size(Y,1));

tSmall=min(per1,per2);
tBig=max(per1,per2);

for ii=1:size(Y,1)
    switch fitmethod
        case 1
            [fit_result,~]=nonlinfit(zts,Y(ii,:) ,ii_best,jj_best);
            per1loc=fit_result.per1; 
            per2loc=fit_result.per2;
        case 2
            beta=nlinfit(zts,Y(ii,:),modelfun,[0 0 0 0 ii_best jj_best]);
            per1loc=beta(5);
            per2loc=beta(6);
        case 3
            beta=fitnlm(zts',Y(ii,:),modelfun,[0 0 0 0 ii_best jj_best]);
            per1loc=table2array(beta.Coefficients('beta5','Estimate'));
            per2loc=table2array(beta.Coefficients('beta6','Estimate'));
    end
    
    % make per1 the smaller one
    
    perSmall_vals(ii)=min(per1loc,per2loc);
    perBig_vals(ii)=max(per1loc,per2loc);
end
period_error=0.5*abs((perSmall_vals-tSmall)/tSmall) + abs((perBig_vals-tBig)/tBig);

end