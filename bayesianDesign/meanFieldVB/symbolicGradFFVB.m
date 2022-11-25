%% calculate fully symbolic gradient
clear
%%
syms S1 S2 S3 
NgaussA=2;NgaussB=2;NgaussT=3;
muvecA=sym('muA',[1 NgaussA]);
muvecB=sym('muB',[1 NgaussB]);
muvecT=sym('muT',[1 NgaussT]);
sdvecA=sym('sdA',[1 NgaussA]);
sdvecB=sym('sdB',[1 NgaussB]);
sdvecT=sym('sdT',[1 NgaussT]);
wvecA=sym('wA',[1 NgaussA]);
wvecB=sym('wB',[1 NgaussB]);
wvecT=sym('wT',[1 NgaussT]);
assume([muvecA muvecB muvecT ...
        sdvecA sdvecB sdvecT ...
        wvecA  wvecB  wvecT ...
        S1     S2     S3],'real')

sqrw=@(w) w.^2/sum(w.^2);
q=( log((  exp(-(S1-muvecA).^2/2./sdvecA.^2)./(sqrt(2*pi)*abs(sdvecA))  )*sqrw(wvecA)') + ...
log((  exp(-(S2-muvecB).^2/2./sdvecB.^2)./(sqrt(2*pi)*abs(sdvecB))  )*sqrw(wvecB)') + ...
log((  exp(-(S3-muvecT).^2/2./sdvecT.^2)./(sqrt(2*pi)*abs(sdvecT))  )*sqrw(wvecT)') );
simplify(gradient(q,wvecA))
%%
symgradq=matlabFunction(gradient(q,[muvecA, sdvecA, wvecA, ...
            muvecB, sdvecB, wvecB, ...
            muvecT, sdvecT, wvecT])','Vars',[...
            muvecA muvecB muvecT ...
            sdvecA sdvecB sdvecT ...
            wvecA  wvecB  wvecT ...
            S1     S2     S3]);

%%
muvecA=rand(1,NgaussA); muvecAc = num2cell(muvecA);
muvecB=rand(1,NgaussB); muvecBc = num2cell(muvecB);
muvecT=rand(1,NgaussT); muvecTc = num2cell(muvecT);
sdvecA=rand(1,NgaussA); sdvecAc = num2cell(sdvecA);
sdvecB=rand(1,NgaussB); sdvecBc = num2cell(sdvecB);
sdvecT=rand(1,NgaussT); sdvecTc = num2cell(sdvecT);
wvecA =rand(1,NgaussA); wvecAc  = num2cell(wvecA);
wvecB =rand(1,NgaussB); wvecBc  = num2cell(wvecB);
wvecT =rand(1,NgaussT); wvecTc  = num2cell(wvecT);
S1=rand;S2=rand;S3=rand;
symgradq(muvecAc{:}, muvecBc{:}, muvecTc{:}, ...
            sdvecAc{:}, sdvecBc{:}, sdvecTc{:}, ...
            wvecAc{:},  wvecBc{:},  wvecTc{:}, ...
            S1,     S2,     S3)
logqgrad([S1 S2 S3],muvecA,sdvecA,wvecA,muvecB,sdvecB,wvecB,muvecT,sdvecT,wvecT)