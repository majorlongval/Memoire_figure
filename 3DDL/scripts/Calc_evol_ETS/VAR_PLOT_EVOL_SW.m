%% VAR PLOT evol SW plots the evolution of the SW area over multiple values of rc in a graph
clear all; close all;clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
R = 1;
r = 0.2;
alpha1 = (pi/2)-0.1;
alpha2 = pi/3;

phic1 = 0;
phic2 = pi/6;

rc = linspace(0, 2*r/3,11);
for i =1:length(rc)
    area1(i) = calc_area_SW_new(alpha1,R,r,rc(i),phic1);
    area2(i) = calc_area_SW_new(alpha2,R,r,rc(i),phic2);
end
plot(rc,area1,'*r');
hold on;
plot(rc,area2,'*b');
plot(rc,area1,'-r');
plot(rc,area2,'-b');

xlabel('x');
ylabel('y');
grid on;
lh = legend('legend111','legend222');
lh.FontSize = 16;
ax = gca;
ax.XTick = linspace(0,2*r/3,11);
ax.XTickLabel = {'s0','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10'};
saveas(gcf,'compar_area','svg')