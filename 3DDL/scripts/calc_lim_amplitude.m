clear all; close all; clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

% Geometric parameters
alpha = pi/6;
R = 0.7033; 
r = 0.2;
rc = 0.06;
phic = 0;
hc = 0;

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
                calc_coeff(alpha, R, r, rc, phic, hc);
            
% Torseur
m = 2;
g = 9.81;
tx = 0;
ty = 0;
tz = 0;
Mx = 0;
My = 0;
Mz = 0;
t = [tx;ty;tz+m*g;Mx;My;Mz];

% The trajectory
pc = [0; 0; 12];
rx = 3.7417;
rz = 0.7654;
ry = 2.8284;
Kx = rx/ry;
Kz = rz/ry;
phix = 1.1832;
phiy = 0;
phiz = 2.4569;
k2 = [Kx*sin(phix);sin(phiy);Kz*sin(phiz)];
k1 = [Kx*cos(phix);cos(phiy);Kz*cos(phiz)];

omegan = sqrt((tz+m*g)/(pc(3)*m));

% discretising omega
omega = linspace(0,4*omegan,1000);

for i =1:length(ay)
   Lambda    = [0 bx(i) cx(i);ay(i) 0 cy(i);az(i) bz(i) 0];
   lambda    = [ak(i);bk(i);ck(i)];
   N         = [Lambda,[dx(i) 0 0;0 ey(i) 0;0 0 hz(i)]];
   n         = [lambda;dk(i);ek(i);hk(i)];
   v         = (pc'*N+n')*t;
   for j = 1:length(omega)
      u = m*omega(j)^2*(pc'*Lambda+lambda')+t'*N';
      theta = sqrt((u*k1)^2+(u*k2)^2);
      cond(j,i) = v/theta;
   end
end

figure; 
plot(cond(:,1)/R,omega/omegan,'-b');
hold on;
plot(cond(:,2)/R,omega/omegan,'--b');
plot(cond(:,3)/R,omega/omegan,'-r');
plot(cond(:,4)/R,omega/omegan,'--r');
plot(cond(:,5)/R,omega/omegan,'-k');
plot(cond(:,6)/R,omega/omegan,'--k');
axis([0 3 0 3]);
grid on;
xlabel('XLABEL'); ylabel('YLABEL');
legend('legend1','legend2','legend3','legend4','legend5','legend6');
saveas(gca,'example_lim_amp','svg');