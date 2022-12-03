%% calculate fully symbolic gradient
clear
%%
syms A1 B1 T1 
Ngauss=5; % number of Gaussians in approximation
d=3;      % number of paramters in model
muvec=sym('mu',[d Ngauss]);
Lvec=sym('L',[d*(d+1)/2 Ngauss]);
Lmat=cell(1,Ngauss);
for ii=1:Ngauss
    Lmat{ii}=INVvech(Lvec(:,ii))
end

theta=[A1;B1;T1];
ii=1
exp(-(theta-muvec(:,ii))'*inv(Lmat{ii}*Lmat{ii}')*(theta-muvec(:,ii)) )
%%
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


function Avech=vech(A)
d=size(A,1);
Avech=[];
for ii=1:size(A,2)
    for jj=ii:d
        Avech(end+1)=A(jj,ii);
    end
end
end


function Lmat = INVvech(Lvec)
d= (8*length(Lvec) + 1)^(1/2)/2 - 1/2;
Lmat = sym(zeros(d,d));
i0=1;
for i=1:d
    Lmat(i:d,i)=Lvec(i0:i0+d-i);
    i0=i0+d+1-i;
end
end
