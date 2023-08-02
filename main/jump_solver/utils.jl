using Distributions
using Optim
function getPower(t,acro,Amp,freq)
  lambda = Amp^2*cos.(2pi*freq*t.-acro)'*cos.(2pi*freq*t.-acro);
  N      = length(t);
  alpha  = 0.05;
  f0=quantile(FDist(2,N-3),1-alpha)
  return 1-cdf(NoncentralF(2,N-3,lambda),f0);
end

function getMinPower(t,Amp,freq)
  power(phi)=getPower(t,phi,Amp,freq);
  op=optimize( power, 0,2pi);
  return op.minimum
end

function get_color_now(ii,N,col1,col2)
# aesthetic for nicer plotting
  rgb= col1 + (ii-1)*(col2-col1)/(N-1)
  return RGB(rgb[1],rgb[2],rgb[3])
end
