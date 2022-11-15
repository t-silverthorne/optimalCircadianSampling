syms ZA Zphi Zf
syms A phi f t S1 S2
syms xiA xiphi xif
syms muA muphi muf
syms sigA sigphi sigf

% priors
piA= exp(-(A-muA)/2/sigA^2)/ZA;
piphi= exp(-(phi-muphi)/2/sigphi^2)/Zphi;
pif= exp(-(f-muf)/2/sigf^2)/Zf;

% get new mean and sdev for q(A)
c=coeffs(collect((A-muA)/2/sigA^2 + (A^2*S1 -2*A*S2),A),A); % no negative
c=fliplr(c);
h=-c(2)/2/c(1);
k=c(3)-c(2)^2/(4*c(1));
c(1)*(A-h)^2+k; % same ply

sd_A=simplify(1/sqrt(2*c(1)));
mu_A=h;

piA_new = exp(-(A-mu_A)/2/sd_A^2)/ZA

%%
%+ log(exp(-(A^2*S1 -2*A*S2) )) )

