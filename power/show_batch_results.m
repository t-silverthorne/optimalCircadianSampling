clear
close all
load('results/Numfreq_10_min_1_max_10_Numamp_8_min_1_max_4.mat')
[X,Y]=meshgrid(ampvals,freqvals);

cmap=[0 0.2 .13].*linspace(1,0,100)'+[0.94 .88 .19].*linspace(0,1,100)'
colormap(cmap)
[M,c]=contourf(X,Y,nu_mat-unif_mat,10,'LineColor','none')
c.LineWidth=2;
