using Base: namemap
using Plots, Plots.Measures 
using LaTeXStrings
include("utils.jl")
#%# random plot example
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2, 
    framestyle=:box, 
    label=nothing, 
    grid=false
)

x=1:0.5:20;
y=1:0.5:10;
f(x,y)=begin
  (3x+y^2)*abs(sin(x)+cos(y));
end
X=repeat(reshape(x,1,:),length(y),1);
Y=repeat(y,1,length(x));
Z=map(f,X,Y);
levels=0:10:200
p1=contourf(x,y,f,fill=true,color=:turbo,cbar=:right,
        titlefontsize=12,
        guidefontsize=12,
        tickfontsize=12
)
plot!(size=(800,600))
savefig("test.pdf")
#%# Actual plot
N=8
namp  = 10;
nfreq = 10;

t     = LinRange(0,1,N+1)
t     = t[1:end-1]
lamp  = LinRange(-2,2,namp);
amp   = [10^i for i in lamp];
freq  = LinRange(1,16,nfreq) 
f(amp,freq) = getMinPower(t,amp,freq)

contourf(amp,freq,f,fill=true,color=:turbo,linewidth=0)
