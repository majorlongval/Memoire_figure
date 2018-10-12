function [pos,vit, acc, t] = traj3Dtrans(pc, a,b,theta, phi, psi, omega, alpha, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%  Defining time
tau = linspace(0,1,1000);
t1 = tau*T;
t2 = linspace(0,(6*pi),1000); 
t3 = t1;
t  = [t1,t1(end)+t2,(t1(end)+t2(end)+t3)];

%% Defining the rotation matricies

% Around z
c1 =cos(theta);
s1 = sin(theta);
Qz = [c1 s1 0;
      -s1 c1 0;
      0    0 1];
  
% Around y
c2 = cos(phi);
s2 = sin(phi);
Qy = [c2  0 s2;
      0   1 0;
      -s2 0 c2];

% Around x
c3 = cos(psi);
s3 = sin(psi);
Qx = [1  0 0;
      0 c3 -s3;
      0 s3 c3];
  
  
%% The equations derived
Qt = Qx*Qy*Qz;
m1 = Qt(1,1);
m2 = Qt(1,2);
m4 = Qt(2,1);
m5 = Qt(2,2);
m7 = Qt(3,1);
m8 = Qt(3,2);

px = atan2(m1*a*cos(alpha)+m2*b*sin(alpha),m2*b*cos(alpha)-m1*a*sin(alpha));
py = atan2(m4*a*cos(alpha)+m5*b*sin(alpha),m5*b*cos(alpha)-m4*a*sin(alpha));
pz = atan2(m7*a*cos(alpha)+m8*b*sin(alpha),m8*b*cos(alpha)-m7*a*sin(alpha));

A11 = m2*b*cos(alpha)-m1*a*sin(alpha);
A12 = m1*a*cos(alpha)+m2*b*sin(alpha);
A21 = m5*b*cos(alpha)-m4*a*sin(alpha);
A22 = m4*a*cos(alpha)+m5*b*sin(alpha);
A31 = m8*b*cos(alpha)-m7*a*sin(alpha);
A32 = m7*a*cos(alpha)+m8*b*sin(alpha);

rx = sqrt(A11^2+A12^2);
ry = sqrt(A21^2+A22^2);
rz = sqrt(A31^2+A32^2);

% setting the initial transition
delta1 = (6*tau.^5-15*tau.^4+10*tau.^3);
delta2 = (30*tau.^4-60*tau.^3+30*tau.^2)/T;
delta3 = (120*tau.^3-180*tau.^2+60*tau)/(T^2);
delta1y = ry*delta1;
delta1x = rx*delta1;
delta1z = rz*delta1;
delta2x = rx*delta2;
delta2y = ry*delta2;
delta2z = rz*delta2;
delta3x = delta3*rx;
delta3y = delta3*ry;
delta3z = delta3*rz;
% setting the full trajectory
delta1y = [delta1y,ry*ones(1,length(t2))];
delta1x = [delta1x,rx*ones(1,length(t2))];
delta1z = [delta1z,rz*ones(1,length(t2))];
delta2x = [delta2x,zeros(1,length(t2))];
delta2y = [delta2y,zeros(1,length(t2))];
delta2z = [delta2z,zeros(1,length(t2))];
delta3x = [delta3x,zeros(1,length(t2))];
delta3y = [delta3y,zeros(1,length(t2))];
delta3z = [delta3z,zeros(1,length(t2))];
% setting the inverse transition trajectory
delta1y = [delta1y,ry*fliplr(delta1)];
delta1x = [delta1x,rx*fliplr(delta1)];
delta1z = [delta1z,rz*fliplr(delta1)];
delta2x = [delta2x,rx*fliplr(delta2)];
delta2y = [delta2y,ry*fliplr(delta2)];
delta2z = [delta2z,rz*fliplr(delta2)];
delta3x = [delta3x,rx*fliplr(delta3)];
delta3y = [delta3y,ry*fliplr(delta3)];
delta3z = [delta3z,rz*fliplr(delta3)];

pos(1:3,:)   = pc+ [delta1x.*sin(omega*t+px);...
    delta1y.*sin(omega*t+py);delta1z.*sin(omega*t+pz)];
vit(1:3,:)  = [delta2x.*sin(omega*t+px)+(omega.*delta1x.*cos(omega*t+px));
               delta2y.*sin(omega*t+py)+(omega.*delta1y.*cos(omega*t+py));
               delta2z.*sin(omega*t+pz)+(omega.*delta1z.*cos(omega*t+pz))];
acc(1:3,:) = [(delta3x-delta1x.*omega^2).*sin(omega*t+px)+...
    2*delta2x.*omega.*cos(omega*t+px);...
    (delta3y-delta1y.*omega^2).*sin(omega*t+py)+...
    2*delta2y.*omega.*cos(omega*t+py);
    (delta3z-delta1z.*omega^2).*sin(omega*t+py)+...
    2*delta2z.*omega.*cos(omega*t+pz)];
end

