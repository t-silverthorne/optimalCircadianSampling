param.useGPU=false;
switch simtype
    case 'rough'
        param.Nperm=1e2;
        param.Nresidual=30; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
        param.maxIter=100;
    case 'fast'
        param.Nperm=1e2;
        param.Nresidual=1e2; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
        param.maxIter=200;
    case 'medium' 
        param.Nperm=1e2;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
        param.maxIter=700;
    case 'long'
        param.Nperm=1e3;
        param.Nresidual=1e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
        param.maxIter=1e3;
    case 'verylong'
        param.Nperm=1e3;
        param.Nresidual=5e3; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples
        param.maxIter=2e3;
end
