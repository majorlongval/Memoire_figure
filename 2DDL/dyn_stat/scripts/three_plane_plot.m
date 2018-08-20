% 3 plane generators
clear all; close all; clc;
addpath('functions');
x  = 5;
y  = 0;
l  = 1;
L  = 5;
ay = 0;
cx = 0;
cy = 0;

geo = pack_geo(ay,cy,cx,L,l);

[t1, t2, t3, t4] = Effort_poly(x, y, geo);
set(gca,'TickLabelInterpreter','Latex')

title('004','FontWeight','normal');
% text(-0.2,0.2,-0.05,'005');
% text(0.4,0,-0.05,'006');
% text(0.4,0,-0.4,'007');
print('torseur_poly','-dsvg');