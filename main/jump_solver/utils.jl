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

function diffMinLambdaDt(t,freq)
# get of min eig wrt measurement schedule
  x1 = cos.(2*pi*freq*t);
  x2 = sin.(2*pi*freq*t);
  X  = [x1 x2]

  D=eigvals(X'*X);
  W=eigvecs(X'*X);
  emin=findmin(D)[2];  # index of min eigenvalue
  v = W[:,emin];
  v = v/norm(v)

  t3  = reshape(t,1,1,length(t));
  pf2 = 2*pi*freq;
  dX  = pf2*[2*cos.(pf2*t3).*sin.(pf2*t3)      cos.(pf2*t3).^2-sin.(pf2*t3).^2;
         cos.(pf2*t3).^2-sin.(pf2*t3).^2  -2*cos.(pf2*t3).*sin.(pf2*t3)];
  return [v'*dX[:,:,ii]*v for ii=1:length(t)]
end
