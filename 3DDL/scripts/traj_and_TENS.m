 clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

R = 0.7033;
r = 0.2;
alpha = pi/6;
g = 9.81;

phic = 0;
hc = 0;
m  = 0.106+0.174+0.523;
rc = 0.1*0.2/m;
tx = 0;
Qx = [1 0 0;0 cos(tx) -sin(tx);0 sin(tx) cos(tx)];
ty = pi/12;
Qy = [cos(ty) 0 sin(ty);0 1 0;-sin(ty) 0 cos(ty)];
tz = pi/6;
Qz = [cos(tz) -sin(tz) 0;sin(tz) cos(tz) 0;0 0 1];

a   = 0.75;
b   = 0.75;
phi = 0;

Q = (Qx*Qy*Qz);

[rx,ry,rz,phix,phiy,phiz] = abQ2rphi(a,b,Q,phi);

pc = [0;0;2];

omega = sqrt(g/pc(3));
T = 20;
nbt = 3;
rat = 0.5;
[t,pos,acc,it1,it2,it3] = calc_ell_traj_full(pc,omega,-rat*rx,-rat*ry,0,phix,phiy,phiz,T,nbt);
% figure;
% plot(t,pos(1,:),t,pos(2,:),t,pos(3,:)); hold on;
plot3(pos(1,:),pos(2,:),pos(3,:));
hold on;
plot3(pos(1,1),pos(2,1),pos(3,1),'*r');
axis([-2 2 -2 2 pc(3)-1 pc(3)+1]);
xlabel('xlabel'); ylabel('ylabel'); zlabel('zlabel'); grid on; hold on;


T = [0;0;0;0;0;0];
tens = calc_tension_traj3D(m,R,r,alpha,rc,phic,hc,T,pos,acc);
figure;
plot(t,tens); grid on;
