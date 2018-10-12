%% Script to explore the effect of the parameters on  the SW with torseur.
clear all; close all; clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
R = 1;
r = 0.2;
alpha = linspace(0,pi/2,25);
rc = 0.05;
pc = 0;
hc = 0;
m  = 1;
tx = 0;
ty = 0;
tz = 0;
Mx = linspace(0,5,25);
My = 0;
Mz = 0;
z_max = 10;
for j =1:length(Mx)
    for i =1:length(alpha)
        vol(i,j) = calc_vol_ETST(alpha(i),R,r,rc,pc,hc,tx,ty,tz,m,Mx(j),My,Mz,z_max);
        fprintf('%d',i);
    end
    fprintf('%d',j);
end

plot(alpha,vol);