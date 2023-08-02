# For comparing Julia code with Matlab code
#
include("utils.jl")
using LinearAlgebra
Nmeas=12;
Amp  =2.23;
freq =4.8

t =[0:.1:1;]

# code for making reduced X

function diffMinLambdaDt(t,freq)
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
