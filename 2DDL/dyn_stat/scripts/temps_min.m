clear all; close all; clc;
L = 5;
l = 0.8;
ay = 0.2;
cy = -0.2;
cx = 0.1;
x0 = 8;
omega = 1.1;
y0 = 2;
rx = 5.48;
ry = 9.75;
phix = -1.08;
phiy = 0.13;
g = 9.81; %m/s²
tors(1) = 0;
tors(2) = 0;
tors(3) = 0;
m =1;
k = rx/ry;
omegan = sqrt(g/x0);
A_1 = -L;
A_2 = ay*(cy-L)+(2*L*cy)-l*(ay+L);
A_3 = -ay*(cy-L)-(2*L*cy)-l*(ay+L);
B_1 = 1;
B_2 = ay-l;
B_3 = -(ay+l);
C_1 = 0;
C_2 = cx*(2*L + ay);
C_3 = -C_2;
D_1  = 0;
D_2  = (2*L+ay);
D_3  = -D_2;
Psi_1 = A_1 + B_1*y0;
Psi_2 = A_2 + B_2*y0;
Psi_3 = A_3 + B_3*y0;
Omega_1 = C_1 + B_1*x0;
Omega_2 = C_2 + B_2*x0;
Omega_3 = C_3 + B_3*x0;
Phi_1 = m*(omega^2)*Psi_1-tors(2)*B_1;
Phi_2 = m*(omega^2)*Psi_2-tors(2)*B_2;
Phi_3 = m*(omega^2)*Psi_3-tors(2)*B_3;
Upsi_1 = ((m*g+tors(1))*B_1-m*omega^2*Omega_1);
Upsi_2 = ((m*g+tors(1))*B_2-m*omega^2*Omega_2);
Upsi_3 = ((m*g+tors(1))*B_3-m*omega^2*Omega_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calcul pour rlimx
d_11 = k*Phi_1*cos(phix) + Upsi_1*cos(phiy);
d_12 = k*Phi_2*cos(phix) + Upsi_2*cos(phiy);
d_13 = k*Phi_3*cos(phix) + Upsi_3*cos(phiy);

d_21 = k*Phi_1*sin(phix) + Upsi_1*sin(phiy);
d_22 = k*Phi_2*sin(phix) + Upsi_2*sin(phiy);
d_23 = k*Phi_3*sin(phix) + Upsi_3*sin(phiy);

theta_1 = sqrt(d_11^2+d_21^2);
theta_2 = sqrt(d_12^2+d_22^2);
theta_3 = sqrt(d_13^2+d_23^2);

A1 = (ry*theta_1)+g*Psi_1;
A2 = (ry*theta_2)-g*Psi_2;
A3 = (ry*theta_3)-g*Psi_3;
% 
B1 = 15*ry*(abs(Omega_1)+abs(Psi_1)+2*ry*k*abs(B_1))/8;
B2 = 15*ry*(abs(Omega_2)+abs(Psi_2)+2*ry*k*abs(B_2))/8;
B3 = 15*ry*(abs(Omega_3)+abs(Psi_3)+2*ry*k*abs(B_3))/8;

C1 = 10*sqrt(3)*ry*(abs(Omega_1)+abs(Psi_1));
C2 = 10*sqrt(3)*ry*(abs(Omega_2)+abs(Psi_2));
C3 = 10*sqrt(3)*ry*(abs(Omega_3)+abs(Psi_3));

T1 = roots([A1 B1 C1]);
T2 = roots([A2 B2 C2]);
T3 = roots([A3 B3 C3]);