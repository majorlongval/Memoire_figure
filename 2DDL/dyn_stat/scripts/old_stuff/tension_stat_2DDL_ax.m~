%% Script pour tracer l'espace de travail statique lorsque ax !=0
clear all;close all;clc;
addpath('functions');



% Setting the geometric parameters of the robot

L = 5;
l = 1;
cy = 0;
cx = 0;
ay = 0;
ax = -0.1;

gp = struct('L',L,'l',l,'cy',cy,'ay',ay,'cx',cx,'ax',ax);

a1 = [gp.ax;gp.ay];
a2 = [0;-gp.l];
a3 = [0;gp.l];
b1 = [0;-gp.L];
b2 = [0;gp.L-gp.l];
b3 = [0;gp.L+gp.l];
c  = [gp.cx;gp.cy];


g = 9.81;

E = [0,-1;1,0];
gam = [g;0;0];

% Creating the surface
x = linspace(0.01,2*L,500);
y = linspace(-(L+l),L,500);


% Calculating the tensions in the cables
for i =1:length(x)
    for j=1:length(y)
        pos(:,j,i) = [x(i);y(j)];
        tau(:,j,i) = tension_stat(pos(:,j,i),gp);
    end
end
tau1(:,:) = tau(1,:,:);
tau2(:,:) = tau(2,:,:);
tau3(:,:) = tau(3,:,:);


% 
% for i =1:length(x)
%     ystat(i,:) = dom_stat([x(i),y(3)],gp);
% end
% figure
% for i =1:length(x)
%     ystatmin(i) = 100;
%     for j =1:length(ystat(1,:))
%         if abs(ystat(i,j))<abs(ystatmin(i))
%             ystatmin(i)=ystat(i,j);
%         end
%     end
% end
% 
% 
% subplot(2,1,1)
% plot(x,ystatmin,'r');
% hold on;
% camroll(-90);
% grid minor;
% plot([x(itest),x(itest)],[-6,6],':k');
% plot([0,5],[ystatmin(itest),ystatmin(itest)],':k');
% plot([0,5],[L,L],'-b');
% axis([0 10 -6 6]);
% subplot(2,1,2)
% hold on
% plot(y,tau(1,:,itest)/g);
% plot(y,tau(2,:,itest)/g);
% plot(y,tau(3,:,itest)/g);
% grid minor
% plot([-6,6],[0,0],':k');
% plot([ystatmin(itest),ystatmin(itest)],[-4,12],':k');
% legend('c1','c2','c3');
% axis([-6,6,-2,2]);
% 
