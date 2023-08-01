using Plots, Plots.Measures 
using LaTeXStrings
include("utils.jl")
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2, 
    framestyle=:box, 
    label=nothing, 
    grid=false,
    titlefontsize=5,
    guidefontsize=5,
    tickfontsize=5,
    legendfontsize=5,
)
##################################################
#%# Fig 2 top: heatmap for uniform and non-uniform
N=8
namp  = 2^6;
nfreq = 2^6;

t     = LinRange(0,1,N+1)
t     = t[1:end-1]
lamp  = LinRange(-2,log10(3),namp);
amp   = [10^i for i in lamp];
freq  = LinRange(1,16,nfreq) 
f(amp,freq) = getMinPower(t,amp,freq)


# left panel
plt1=heatmap(amp,freq,f,clim=(0,1),levels=10,fill=true,color=:viridis,linewidth=0,
       xlabel=L"$A$", ylabel=L"$f$")
plot!(colorbar=false)


# right panel
t=randn(length(t))
plt2=heatmap(amp,freq,f,clim=(0,1),levels=10,fill=true,color=:viridis,linewidth=0,
       xlabel=L"$A$")
plot!(yticks=false)
plot!(colorbar=false,location=:none)


# dummy plot for adding colorbar
h2 = scatter([NaN NaN], [NaN NaN], zcolor=NaN, clims=(0,1),
                 xlims=(0,0), xshowaxis=false, 
                 yshowaxis=false, label="", c=:viridis,
                 colorbar_title=L"\min_\phi \gamma(\phi; A,f)", grid=false,framestyle=:white,background_color=:white,
                 framecolor=:white,foreground_color_border=:white,
                 foreground_color_axis=:white,colorbar_titlefontsize=5)

# combine
layout = @layout [grid(1,2) b{0.03w}]
plt=plot(plt1,plt2,h2,layout=layout)

# export
plot!(size=(480,200))
savefig("~/research/overleaf/rate_limited_sampling/figures/fig1_top.pdf")


##################################################
#%# Fig 2 bottom, line graph for uniform and optimal 
freq_list = LinRange(3,4,nfreq);
t     = LinRange(0,1,N+1)
t     = t[1:end-1]
nfreq=16
pv= [0:.01:2*pi;]

# aesthetic
col1=[1,0,0];
col2=[0,1,0];


# left panel
plt1=plot()
for ii in [1:length(freq_list);]
  fr=freq_list[ii]
  f(phi) = getPower(t,phi,2,fr)
  fv=map(f,pv)
  plt1=plot!(pv,fv,lc=get_color_now(ii,length(freq_list),col1,col2))
  
  display(plt1)
end
plot!(ylims=(0,1),framecolor=:white,framestyle=:white)
plot!(xlabel=L"$\phi$")
plot!(ylabel=L"$\gamma(\phi;A,f)$")

# right panel
t=randn(length(t))
plt2=plot()
for ii in [1:length(freq_list);]
  fr=freq_list[ii]
  f(phi) = getPower(t,phi,2,fr)
  fv=map(f,pv)
  plt2=plot!(pv,fv,lc=get_color_now(ii,length(freq_list),col1,col2))
  
  display(plt2)
end
plot!(ylims=(0,1),framecolor=:white,framestyle=:white)
plot!(xlabel=L"$\phi$")
plot!(yticks=false)


# dummy plot for colorbar
h2=scatter([NaN NaN],[NaN NaN],zcolor=NaN,
           xshowaxis=false, yshowaxis=false,
           clims=(freq_list[1],freq_list[end]),
           c=cgrad([RGB(col1[1],col1[2],col1[3]),
           RGB(col2[1],col2[2],col2[3])],
          [0,1]),colorbar=true,
      location=:left,colorbar_titlefontsize=5,
      colorbar_title=L"$f$")

# combine
layout = @layout [grid(1,2) b{0.03w}]
plt=plot(plt1,plt2,h2,layout=layout)

# export
plot!(size=(480,200))
savefig("~/research/overleaf/rate_limited_sampling/figures/fig1_bottom.pdf")
