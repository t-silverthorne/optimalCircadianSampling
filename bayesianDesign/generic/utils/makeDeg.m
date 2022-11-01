function pvec = makeDeg(tilde_pvec,model)
switch model
    case 'cosinorOneFreq'
        pvec(1)=tilde_pvec(1); % positive amplitude
        pvec(2)=atanh(tilde_pvec(2)/pi-1); % acro in [0,2pi]
        pvec(3)=tilde_pvec(3); % period in [0,1]
end