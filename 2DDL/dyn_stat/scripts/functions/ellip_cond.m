function [omega, omegan, rlimx_1, rlimx_2, rlimx_3,...
    rlimy_1, rlimy_2, rlimy_3] = ...
    ellip_cond(L, l, cx, cy, ay, x0, y0,phix,phiy, k, tors,m,res)
%ELLIP_COND takes the geo. parameters of the planar mechanism and the wrench and shoots
%out a graphic showing the surfaces of accesible range of motion.
%   The first 7 parameters are well know. The last one is the ratio between
%   rx and ry such that ry = krx.
g = 9.81; %m/s²
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_1 = A_1 + B_1*y0;
Psi_2 = A_2 + B_2*y0;
Psi_3 = A_3 + B_3*y0;

Omega_1 = C_1 + B_1*x0;
Omega_2 = C_2 + B_2*x0;
Omega_3 = C_3 + B_3*x0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    omega = linspace(0,4*omegan,res);
for i =1:length(omega)    
    Phi_1 = m*(omega(i)^2)*Psi_1-tors(2)*B_1;
    Phi_2 = m*(omega(i)^2)*Psi_2-tors(2)*B_2;
    Phi_3 = m*(omega(i)^2)*Psi_3-tors(2)*B_3;
    
    Upsi_1 = ((m*g+tors(1))*B_1-m*omega(i)^2*Omega_1);
    Upsi_2 = ((m*g+tors(1))*B_2-m*omega(i)^2*Omega_2);
    Upsi_3 = ((m*g+tors(1))*B_3-m*omega(i)^2*Omega_3);
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rlimy_1(i) = -(Psi_1*(m*g+tors(1))-Omega_1*tors(2)+D_1*tors(3))*...
                   theta_1/(d_11^2+d_21^2);
    rlimy_2(i) = -(Psi_2*(m*g+tors(1))-Omega_2*tors(2)+D_2*tors(3))*...
                   theta_2/(d_12^2+d_22^2);
    rlimy_3(i) = -(Psi_3*(m*g+tors(1))-Omega_3*tors(2)+D_3*tors(3))*...
                   theta_3/(d_13^2+d_23^2);
    
    rlimx_1(i) = k*rlimy_1(i);
    rlimx_2(i) = k*rlimy_2(i);
    rlimx_3(i) = k*rlimy_3(i);
end

end
