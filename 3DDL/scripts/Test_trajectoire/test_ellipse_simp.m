%% Ce script à pour but de tester la forme simplifiée des équations
clear all; close all; clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
g = 9.8066;

pc = [0;0;12];
a  =    2;
b  =    2;
theta = pi/4;
phi =   pi/4;
psi =   0;
omega = sqrt(g/pc(3));
decal = 0;
[pos,vit,acc,t]=traj3Dtrans(pc, a,b,theta, phi, psi, omega, decal, 20*pi/omega)
%comet3(pos(1,:),pos(2,:),pos(3,:));

R = 0.7033;
r = 0.2;
g = 9.8066;
m = 1;
c = [0.06;0;0];
tx = 0;
ty = 0;
tz = 0;
Mx = 0;
My = 0;
Mz = 0;
a11 = r*[0; 1; 0];
a12 = -a11;
ct = cos(2*pi/3); st = sin(2*pi/3);
Qs = [ct -st 0;st ct 0;0 0 1];
a21 = Qs*a11;
a22 = -a21;
a31 = Qs*a21;
a32 = -a31;
alpha = pi/6;
R1 = R*[cos(alpha);sin(alpha);0];
R2 = Qs*R1;
R3 = Qs*R2;






b11 = R1+a11;
b12 = R1+a12;
b21 = R2+a21;
b22 = R2+a22;
b31 = R3+a31;
b32 = R3+a32;

for i =1:length(pos)
rho1(:,i) = pos(:,i)-R1;
rho2(:,i) = rho1(:,i);
rho3(:,i) = pos(:,i)-R2;
rho4(:,i) = rho3(:,i);
rho5(:,i) = pos(:,i)-R3;
rho6(:,i) = rho5(:,i);
e1(:,i) = rho1(:,i)/norm(rho1(:,i));
e2(:,i) = rho3(:,i)/norm(rho3(:,i));
e3(:,i) = rho5(:,i)/norm(rho5(:,i));
delta11 = cross(a11-c,e1(:,i));
delta12 = cross(a12-c,e1(:,i));
delta21 = cross(a21-c,e2(:,i));
delta22 = cross(a22-c,e2(:,i));
delta31 = cross(a31-c,e3(:,i));
delta32 = cross(a32-c,e3(:,i));

M(:,:,i) = [rho1(:,i)/norm(rho1(:,i)) rho2(:,i)/norm(rho2(:,i)) rho3(:,i)/norm(rho3(:,i))...
            rho4(:,i)/norm(rho4(:,i)) rho5(:,i)/norm(rho5(:,i)) rho6(:,i)/norm(rho6(:,i));
     delta11 delta12 delta21 delta22 delta31 delta32];

gamma(:,i) = [tx-m*acc(1,i);ty-m*acc(2,i);tz-m*acc(3,i)+m*g;Mx;My;Mz];
tens(:,i) = M(:,:,i)\gamma(:,i);
drho(:,i) = (M(:,:,i)')\[vit(:,i);0;0;0];
end


for i =1:6
    for j =1:length(drho(i,:))
        if drho(i,j)>0
            drho_bin(i,j)=1;
        elseif drho(i,j)<0
            drho_bin(i,j)=-1;
        else
            drho_bin(i,j)=0;
        end
            
    end
end



plot(t, tens(1,:));
hold on
plot(t, tens(2,:));
plot(t, tens(3,:));
plot(t, tens(4,:));
plot(t, tens(5,:));
plot(t, tens(6,:));
% 
% 
% torque1  = (tens(1,:)+tens(2,:))*0.0155;
% torque2  = (tens(3,:)+tens(4,:))*0.0155;
% torque3  = (tens(5,:)+tens(6,:))*0.0155;
% 
% figure
% plot(t,torque1);
% hold on;
% plot(t,torque2);
% plot(t,torque3);
% 
% figure
% 
% plot(t, drho_bin(1,:));
% hold on
% plot(t, drho_bin(2,:));
% plot(t, drho_bin(3,:));
% plot(t, drho_bin(4,:));
% plot(t, drho_bin(5,:));
% plot(t, drho_bin(6,:));
% 
% 
% figure
% for i =1:6
% plot(t, 0.02*drho(i,:).*tens(i,:)./abs(drho(i,:)));
% hold on;
% end
% for i =1:6
%     avg_tens(i) = mean(tens(i,:));
% end