function theta=getTheta(tvec,fnames)
tc=num2cell(tvec);
theta=cell2struct(tc(:),fnames);
end