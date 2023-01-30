param.useGPU=false;
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e2; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'medium' 
        param.Nperm=1e2;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'long'
        param.Nperm=1e3;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
    case 'verylong'
        param.Nperm=1e3;
        param.Nresidual=5e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
end
