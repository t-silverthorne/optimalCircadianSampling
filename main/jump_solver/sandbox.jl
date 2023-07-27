#%# Functions
using LinearAlgebra
function constructX(t,p)
# Construct the design matrix 
    x0=ones(length(t),1);
    x1=sin.(2*pi*p.freq*t);
    x2=cos.(2*pi*p.freq*t);
    return [x0 x1 x2];
end

function constructReducedX(t,p)
# Construct the design matrix assuming beta0=0 
    x1=sin.(2*pi*p.freq*t);
    x2=cos.(2*pi*p.freq*t);
    return [x1 x2];
end

Base.@kwdef mutable struct harmonic_param
# store parameters for harmonic regression 
    Amp::Float64=1;
    freq::Float64=1;
    phase::Float64=0;
end


#%# Main program
p=harmonic_param();
t=[0:.1:1;];
t=t[1:end-1];

X=constructX(t,p);

