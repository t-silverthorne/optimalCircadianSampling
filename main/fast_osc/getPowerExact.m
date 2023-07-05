function beta=getPowerExact(Amp,acro,freq,mt)
Nmeas    = length(mt); % num samples
alpha    = 0.05; % type I error rate

csq    = cos(2*pi*freq*mt-acro);
lambda   =  Amp^2*csq*csq'; % non-centrality

beta=1-ncfcdf( finv(1-alpha,2,Nmeas-3) ,2,Nmeas-3,lambda);
end