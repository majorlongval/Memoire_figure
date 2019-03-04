%% Calculating the transition period for elliptical trajectory
clear all; close all; clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L      = 5;            % m distance between independant and dep. pulleys
l      = 0.8;          % m width of end-effector
cx     = 0.2;            % m X position of centre of mass
cy     = -0.1;            % m Y position of centre of mass
ay     = 0.3;          % m attach position of independant cable
x0     = 2.5;           % m X centre of oscillation
y0     = 1;          % m Y centre of oscillation
m      = 2;            % kg mass of the end-effector
mC     = 1;
dx     = 0.4;
dy     = -0.2;
ex     = 0.3;
ey     = -0.3;
res    = 1000;         % resolution of functions
g      = 9.81;% m/s^2 gravitational acceleration
tors   = mC*g*[1, 0 , dy];   % [N, N, Nm] The wrench

omega = 1.86; %rad/s


ry = 8;
k = 1/16;
rx = ry*k;
phix = 0;
phiy = 0;

% Parameters
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
d_11 = k*Phi_1*cos(phix) + Upsi_1*cos(phiy);
d_12 = k*Phi_2*cos(phix) + Upsi_2*cos(phiy);
d_13 = k*Phi_3*cos(phix) + Upsi_3*cos(phiy);

d_21 = k*Phi_1*sin(phix) + Upsi_1*sin(phiy);
d_22 = k*Phi_2*sin(phix) + Upsi_2*sin(phiy);
d_23 = k*Phi_3*sin(phix) + Upsi_3*sin(phiy);

theta_1 = sqrt(d_11^2+d_21^2);
theta_2 = sqrt(d_12^2+d_22^2);
theta_3 = sqrt(d_13^2+d_23^2);

kappa_1 = Psi_1*(m*g+tors(1))-Omega_1*tors(2)+D_1*tors(3);
kappa_2 = Psi_2*(m*g+tors(1))-Omega_2*tors(2)+D_2*tors(3);
kappa_3 = Psi_3*(m*g+tors(1))-Omega_3*tors(2)+D_3*tors(3);


% Calculating the times
A1 = (10*sqrt(3)/3)*(abs(Omega_1)+k*Psi_1);
A2 = (10*sqrt(3)/3)*(abs(Omega_2)+k*Psi_2);
A3 = (10*sqrt(3)/3)*(abs(Omega_3)+k*Psi_3);

B1 = (15/8)*(abs(Omega_1)+k*Psi_1+2*k*abs(B_1));
B2 = (15/8)*(abs(Omega_2)+k*Psi_2+2*k*abs(B_2));
B3 = (15/8)*(abs(Omega_3)+k*Psi_3+2*k*abs(B_3));

C1 = ry*theta_1+kappa_1;
C2 = ry*theta_2+kappa_2;
C3 = ry*theta_3+kappa_3;

T1 = roots([A1 B1 C1]);
T2 = roots([A2 B2 C2]);
T3 = roots([A3 B3 C3]);
Tall = [T1;T2;T3];
T = max(Tall);