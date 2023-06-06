close all
clear all

n=3; % number of columns
m=2; % number of rows
tiledlayout(m,n)
ax1=nexttile(1)
xlim([0 1])
ylim([2 3])

ax2=nexttile(2)
xlim([-1 1])
ylim([1 3])

ax3=nexttile(3)
xlim([.1 .21])
ylim([1 2])

ax4=nexttile(4)
xlim([-5 5])
ylim([1 3])

ax5=nexttile(5)
xlim([0 1])
ylim([1 2])

ax6=nexttile(6)
xlim([-1 1])
ylim([1 1.1])

linkaxes([ax1 ax2 ax3],'y')

linkaxes([ax1 ax4],'x')
