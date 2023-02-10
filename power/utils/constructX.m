function X=constructX(t,param)
if isfield(param,'freq_est')
    freq=param.freq_est;
else
    freq=param.freq_true;
end
if param.useGPU
    x0=gpuArray(ones(1,length(t)));
    x1=gpuArray(sin(2*pi*freq*t));
    x2=gpuArray(cos(2*pi*freq*t));
else
    x0=ones(1,length(t));
    x1=sin(2*pi*freq*t);
    x2=cos(2*pi*freq*t);
end
X= [x0' x1' x2'];

end