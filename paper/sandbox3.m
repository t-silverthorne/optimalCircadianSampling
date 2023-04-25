figure
%%
clf
numsc=4;
c1rgb=[100, 143, 255]/256; % uniform colour
c2="#785EF0"; % 2-uniform colour
c3="#DC267F"; % Unif(0,1) colour
c4rgb=[254, 97, 0]/256; % jittered colour

scale_range=[1e-8 1e-4 1e-2 1e-1]
numsc=length(scale_range)


newx=[-5 -4 -3 -2];
colors=[c1rgb;c4rgb];
oldx=[-5 -2]
cloc=interp1(oldx,colors,newx)

for jj=1:numsc
    t_jitloc = t_unif + scale_range(jj)*rand(1,N);
    [pwr_est,amp_est,phi_est] = getPowerBatch(t_jitloc,p,I3,I4);
    histogram(amp_est(:,acind),'EdgeColor','none','FaceColor',cloc(jj,:),'Normalization','probability')
    hold on
    drawnow
end
colormap(cloc)