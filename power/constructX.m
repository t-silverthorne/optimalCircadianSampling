function X=constructX(t,param)
per=param.per;
if param.useGPU
    x0=gpuArray(ones(1,length(t)));
    x1=gpuArray(sin(2*pi*per*t));
    x2=gpuArray(cos(2*pi*per*t));
else
    x0=ones(1,length(t));
    x1=sin(2*pi*per*t);
    x2=cos(2*pi*per*t);
end
X= [x0' x1' x2'];

end