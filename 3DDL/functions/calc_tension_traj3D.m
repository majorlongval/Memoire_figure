function tens = calc_tension_traj3D(m,R,r,alpha,rc,phic,hc,t,pos,acc)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
g = 9.81;
c = [rc*cos(phic);rc*sin(phic);hc];
tx = t(1);ty = t(2);tz = t(3); Mx = t(4); My = t(5); Mz = t(6);
a11 = r*[0; 1; 0];
a12 = -a11;
ct = cos(2*pi/3); st = sin(2*pi/3);
Qs = [ct -st 0;st ct 0;0 0 1];
a21 = Qs*a11;
a22 = -a21;
a31 = Qs*a21;
a32 = -a31;
R1 = R*[cos(alpha);sin(alpha);0];
R2 = Qs*R1;
R3 = Qs*R2;

% 
% b11 = R1+a11;
% b12 = R1+a12;
% b21 = R2+a21;
% b22 = R2+a22;
% b31 = R3+a31;
% b32 = R3+a32;

for i =1:length(pos(1,:))
rho1 = pos(:,i)-R1;
rho2 = rho1;
rho3 = pos(:,i)-R2;
rho4 = rho3;
rho5 = pos(:,i)-R3;
rho6 = rho5;
e1 = rho1/norm(rho1);
e2 = rho3/norm(rho3);
e3 = rho5/norm(rho5);
delta11 = cross(a11-c,e1);
delta12 = cross(a12-c,e1);
delta21 = cross(a21-c,e2);
delta22 = cross(a22-c,e2);
delta31 = cross(a31-c,e3);
delta32 = cross(a32-c,e3);

M = [e1 e1 e2 e2 e3 e3;...
     delta11 delta12 delta21 delta22 delta31 delta32];

gamma = [tx-m*acc(1,i);ty-m*acc(2,i);tz-m*acc(3,i)+m*g;Mx;My;Mz];
tens(:,i) = M\gamma;
end

