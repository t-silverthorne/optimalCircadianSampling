function tilde_pvec = makeNonDeg(pvec,model)
switch model
    case 'cosinorOneFreq'
        tilde_pvec(1)=abs(pvec(1)); % positive amplitude
        tilde_pvec(2)=pi+pi*tanh(pvec(2)); % acro in [0,2pi]
        tilde_pvec(3)=abs(pvec(3)); % period in [0,1]
end

