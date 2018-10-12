%% Script to plot all the planes of the limits of the SW with torseur
clear all; close all;clc;

R = 1;
r = 0.2;
alpha = pi/4;
r_c = 0;
phi_c = 0;
h_c = 0;
m  = 1;
g  = 9.81;
tx = 0.5;
ty = 0;
tz = 0;
Mx = 0.5;
My = 0;
Mz = 0;
Tz = m*g+tz;

px(1) = Mx - h_c*ty + Tz*(r_c*sin(phi_c) + 2*r*cos(alpha)^2);
px(2) = h_c*ty - Mx - Tz*(r_c*sin(phi_c) - 2*r*cos(alpha)^2);
px(3) = Mx - Tz*(r*cos(alpha)^2 - r_c*sin(phi_c) + (3^(1/2)*r*sin(2*alpha))/2) - h_c*ty;
px(4) = h_c*ty - Tz*(r_c*sin(phi_c) + r*cos(alpha)^2 + (3^(1/2)*r*sin(2*alpha))/2) - Mx;
px(5) = Mx + Tz*(r_c*sin(phi_c) - r*cos(alpha)^2 + (3^(1/2)*r*sin(2*alpha))/2) - h_c*ty;
px(6) = h_c*ty - Tz*(r_c*sin(phi_c) + r*cos(alpha)^2 - (3^(1/2)*r*sin(2*alpha))/2) - Mx;

py(1) = h_c*tx + Tz*(r*sin(2*alpha) - r_c*cos(phi_c));
py(2) = Tz*(r*sin(2*alpha) + r_c*cos(phi_c)) - h_c*tx;
py(3) = h_c*tx - Tz*((r*sin(2*alpha))/2 + r_c*cos(phi_c) - 3^(1/2)*r*cos(alpha)^2);
py(4) = Tz*(r_c*cos(phi_c) - (r*sin(2*alpha))/2 + 3^(1/2)*r*cos(alpha)^2) - h_c*tx;
py(5) = h_c*tx - Tz*((r*sin(2*alpha))/2 + r_c*cos(phi_c) + 3^(1/2)*r*cos(alpha)^2);
py(6) = - Tz*((r*sin(2*alpha))/2 - r_c*cos(phi_c) + 3^(1/2)*r*cos(alpha)^2) - h_c*tx;

pz(1) = Mz - tx*(2*r - 2*r*sin(alpha)^2 + r_c*sin(phi_c)) - ty*(r*sin(2*alpha) - r_c*cos(phi_c));
pz(2) = tx*(2*r*sin(alpha)^2 - 2*r + r_c*sin(phi_c)) - Mz - ty*(r*sin(2*alpha) + r_c*cos(phi_c));
pz(3) = Mz + ty*((r*sin(2*alpha))/2 + r_c*cos(phi_c) - 3^(1/2)*r*cos(alpha)^2) + tx*(r*cos(alpha)^2 - r_c*sin(phi_c) + (3^(1/2)*r*sin(2*alpha))/2);
pz(4) = tx*(r_c*sin(phi_c) + r*cos(alpha)^2 + (3^(1/2)*r*sin(2*alpha))/2) - ty*(r_c*cos(phi_c) - (r*sin(2*alpha))/2 + 3^(1/2)*r*cos(alpha)^2) - Mz;
pz(5) = Mz + ty*((r*sin(2*alpha))/2 + r_c*cos(phi_c) + 3^(1/2)*r*cos(alpha)^2) - tx*(r_c*sin(phi_c) - r*cos(alpha)^2 + (3^(1/2)*r*sin(2*alpha))/2);
pz(6) = ty*((r*sin(2*alpha))/2 - r_c*cos(phi_c) + 3^(1/2)*r*cos(alpha)^2) - Mz + tx*(r_c*sin(phi_c) + r*cos(alpha)^2 - (3^(1/2)*r*sin(2*alpha))/2);

pk(1) = (-2*R*h_c*cos(alpha))*ty+(R*r*cos(alpha) + 2*R*r_c*cos(alpha)*sin(phi_c))*Tz +...
      Mx*2*R*cos(alpha);
pk(2) = (2*R*h_c*cos(alpha))*ty+(R*r*cos(alpha) - 2*R*r_c*cos(alpha)*sin(phi_c))*Tz -...
      Mx*2*R*cos(alpha);
pk(3) = tx*(3^(1/2)*R*h_c*cos(alpha))+ty*(R*h_c*cos(alpha))+...
      Tz*(R*r*cos(alpha) - R*r_c*cos(alpha)*sin(phi_c) - 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c)) -...
      Mx*R*cos(alpha) + 3^(1/2)*R*cos(alpha)*My;
pk(4) = tx*(-3^(1/2)*R*h_c*cos(alpha))-ty*(R*h_c*cos(alpha))+...
      Tz*(R*r*cos(alpha) + R*r_c*cos(alpha)*sin(phi_c) + 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c)) +...
      Mx*R*cos(alpha) - 3^(1/2)*R*cos(alpha)*My;
pk(5) = tx*(-3^(1/2)*R*h_c*cos(alpha))+ty*(R*h_c*cos(alpha))+...
      Tz*(R*r*cos(alpha) - R*r_c*cos(alpha)*sin(phi_c) + 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c)) -...
      Mx*R*cos(alpha) - 3^(1/2)*R*cos(alpha)*My;
pk(6) = tx*(3^(1/2)*R*h_c*cos(alpha))-ty*(R*h_c*cos(alpha))+...
      Tz*(R*r*cos(alpha) + R*r_c*cos(alpha)*sin(phi_c) - 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c)) +...
      Mx*R*cos(alpha) + 3^(1/2)*R*cos(alpha)*My;




[x y] = meshgrid(-1:0.1:1); % Generate x and y data
z1 = -(1/pz(1)).*(px(1).*x + py(1).*y + pk(1));
z2 = -(1/pz(2)).*(px(2).*x + py(2).*y + pk(2));
z3 = -(1/pz(3)).*(px(3).*x + py(3).*y + pk(3));
z4 = -(1/pz(4)).*(px(4).*x + py(4).*y + pk(4));
z5 = -(1/pz(5)).*(px(5).*x + py(5).*y + pk(5));
z6 = -(1/pz(6)).*(px(6).*x + py(6).*y + pk(6));

surf(x,y,z1,'FaceColor',[1 0 0]) %Plot the surface 1
hold on;
surf(x,y,z2) %Plot the surface 1
surf(x,y,z3) %Plot the surface 1
surf(x,y,z4) %Plot the surface 1
surf(x,y,z5) %Plot the surface 1
surf(x,y,z6) %Plot the surface 1

axis([-1 1 -1 1 0 5]);
set(gca,'Zdir','reverse')