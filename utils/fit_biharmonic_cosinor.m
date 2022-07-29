function betas = fit_biharmonic_cosinor(Y,zts,per1,per2)
% inputs:
% Y = time-series data, each columns correspond to measurement times 
%       and rows correspond to independent samples
% zts = measurement times (must be same length as size(Y,1)
% per1 = first period of harmonic regression model
% per2 = second period of harmonic regression model

x1=sin(2*pi*zts/per1);
x2=cos(2*pi*zts/per1);
x0=ones(1,numel(zts));

x3=sin(2*pi*zts/per2);
x4=cos(2*pi*zts/per2);
X= [x0' x1' x2' x3' x4'];

% do linear regression 
betas=(X'*X)\X'*Y'; 

end

