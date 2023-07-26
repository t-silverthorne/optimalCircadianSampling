function betamin=getMinPower(Amp,freq,mt)
[~,betamin]=fminbnd(@(phi) getPowerExact(Amp,phi,freq,mt),0,2*pi);
end