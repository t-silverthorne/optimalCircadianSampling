function X=constructX(t,param)
freq_est=param.freq_est;
if param.useGPU
    x0=gpuArray(ones(1,length(t)));
    x1=gpuArray(sin(2*pi*freq_est*t));
    x2=gpuArray(cos(2*pi*freq_est*t));
else
    x0=ones(1,length(t));
    x1=sin(2*pi*freq_est*t);
    x2=cos(2*pi*freq_est*t);
end
X= [x0' x1' x2'];

end