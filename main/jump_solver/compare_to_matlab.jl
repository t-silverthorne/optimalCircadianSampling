# imports
using Optim: minimizer
using LinearAlgebra
using Optim
using Printf
include("utils.jl")

# parameters 
Nmeas=8;
Amp  =2.23;
freq =Nmeas*0.9;

# measurement grid 
t = LinRange(0,1,Nmeas+1)'
t = Vector(t[1:Nmeas]);
t = t .+ minimum(diff(t))/2; # penalty methods dont like starting on boundary

# optimizer options 
lower   = zeros(length(t));
upper   = ones(length(t));
initial = t;


# Method 1: directly optimize power
f(t) = -getMinPower(t,Amp,freq)
res1=optimize(f,lower,upper,initial);


# Method 2: optimize lambda, (FD derivative)
f(t) = -getMinLambda(t,freq)
res2=optimize(f,lower,upper,initial);
res2_fval=getMinPower(res.minimizer,Amp,freq)

# Method 3: optimize lambda with user-supplied gradient
function grad_wrapper!(storage,x)
  storage[1:length(t)]=diffMinLambdaDt(x,freq)
end
initial
res3=optimize(f,grad_wrapper!,lower,upper,initial);
res3_fval=getMinPower(res3.minimizer,Amp,freq)

#%# Method 4: use Eoptimality
eOpt=getLambdaEoptMethod(freq,Int(1e2),Nmeas)

#%#
# Summary
@printf "\nUniform:              %2.3f\n" getMinPower(t,Amp,freq) 
@printf "LBFGS naive:          %2.3f   %2.3f\n" res1.minimum  res1.time_run  
@printf "LBFGS lambda:         %2.3f   %2.3f\n" res2_fval     res2.time_run  
@printf "LBFGS diff lambda:    %2.3f   %2.3f\n" res3_fval     res3.time_run
