param.useGPU=true;
switch simtype
    case 'rough'
        param.Nperm=5;
        param.Nresidual=5; % SMALL RIGHT NOW
        param.Nacro=2; % num. fourier samples

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
        param.Nperm=1e4;
        param.Nresidual=2e2; % SMALL RIGHT NOW
        param.Nacro=32; % num. fourier samples

end
