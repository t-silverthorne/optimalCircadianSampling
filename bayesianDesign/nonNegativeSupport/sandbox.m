addpath('utils/')
NL=4; % set up time grid
NR=4;
tauA=0;
tauB=1/3;

[mt_unif,mt_nu]=getSamplingSchedules(NL,NR,tauA,tauB);

settings.model='cosinorOneFreq';
M=getBayesianFIMcirc(NL+NR,settings.model); % just a function nothing evaluated yet

for ii=1:5
    mt_unif_long=repmat(mt_nu,1,ii);
    mt_unifc=num2cell(mt_unif_long);
    M=getBayesianFIMcirc((NL+NR)*ii,settings.model);
    log(det(M(0.1,0.1,2.5,mt_unifc{:})))
end