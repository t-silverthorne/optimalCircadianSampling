addpath('../utils')
%%
t=rand(10);
t=t(1:end-1);
p.Amp=1;
p.acro=0;
p.freq=1;
X=constructX(t,p)
det(X'*X)
%%
syms A B C D E F G H I Nm
assume(c>0)
assume(s>0)
assume(Nm>0)
syms N c s c2 s2 cs
XtX=[N c   s;
     c c2 cs;
     s cs s2];
diag(simplify(inv(XtX)))

%%
XtXinv=inv([A B C; D E F; G H I]);
XtXinv=subs(XtXinv,[A B C],[1 c s])
%%
XtXinv=subs(XtXinv,[D E F],[c c^2 s*c]);
%%
XtXinv=subs(XtXinv,[G H I],[s c*s s^2]);