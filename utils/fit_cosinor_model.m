function cosinor_stats=fit_cosinor_model(Y,zts,per)
% Reimplementation of rowCosinor code from the Petronis lab, performs
% harmonic regression on a matrix of independent time series.
% 
% inputs:
% Y = time-series data, each columns correspond to measurement times 
%       and rows correspond to independent samples
% zts = measurement times (must be same length as size(Y,1)
% per = period of harmonic regression model
%
% outputs:
% for each row - acrophase,amplitude,mesor,rsq,Fstatistic,pvalue

x1=sin(2*pi*zts/per);
x2=cos(2*pi*zts/per);
x0=ones(1,size(Y,2));
X= [x0' x1' x2'];

% do linear regression 
betas=(X'*X)\X'*Y'; 
cosinor_stats.acrophases_rad=atan2(betas(2,:),betas(3,:)); % acrophase in radians
cosinor_stats.acrophases=mod(per/2/pi*cosinor_stats.acrophases_rad,per); % acrophase in hours
cosinor_stats.amplitudes=sqrt(betas(2,:).*betas(2,:)+betas(3,:).*betas(3,:));% modified from definition in R code
cosinor_stats.mesor=betas(1,:);

fits=(X*betas)';

% calculate R^2 value
SStot=sum((Y-mean(Y,2)).^2,2);
SSres=sum((fits-Y).^2,2);
cosinor_stats.rsq=1-(SSres./SStot);

% do an F test
SSmod=SStot-SSres;
DFres = size(Y,2) - 3;
DFmod = 2;
MSres = SSres./DFres;
MSmod = SSmod./DFmod;
cosinor_stats.Fstatistic = MSmod./MSres;

% p-value from F-test
get_pvalue = @(y)  1-integral(@(x) fpdf(x,DFmod,DFres),-Inf,y);
cosinor_stats.pval=arrayfun(get_pvalue,cosinor_stats.Fstatistic);

end
