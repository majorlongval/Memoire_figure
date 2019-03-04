%% Script to determine the tension in time of a certain ellipse
clear all;close all; clc;

addpath('functions');

% GEOMETRY
m =2;
g = 9.81;
L = 5;
l = 0.8;
ay = 0.3;
cy = -0.1;
cx = 0.2;
x0 = 2.5;
y0 = 1;
pos_centre = [x0;y0];
dy = -0.2;
% TRAJECTORY
phix = 0;
phiy = 0;
ry = 8;
k = 1/16;
rx = k*ry;
omega = 1.86;
T = 4*pi/omega;



[pos, acc, t] = ...
    traj2Drxry(pos_centre, rx, ry, omega, phix, phiy);
[posup,posdown, accup,accdown,tup,tdown] = ...
    transtraj2D(pos_centre, rx, ry,phix,phiy, omega, T);


geo = pack_geo(ay,cy,cx,L,l);

mC = 1;

Tors = mC*g*[1;0;dy];

postot = [posup,pos,posdown];
% 

acctot = [accup,acc,accdown];
% 
ttot = [tup tup(end)+t tup(end)+t(end)+tdown];

% plot(t,pos(1,:),t,pos(2,:));
% figure
% plot(t,acc(1,:),t,acc(2,:));

for i =1:length(tup)
    tensup(:,i) = tension(m,geo,posup(:,i),accup(:,i),Tors);
end
for i =1:length(t)
    tens(:,i) = tension(m,geo,pos(:,i),acc(:,i),Tors);
end
for i =1:length(tdown)
    tensdown(:,i) = tension(m,geo,posdown(:,i),accdown(:,i),Tors);
end
figure;
plot(tup,tensup(1,:),'--b',tup,tensup(2,:),'--r',tup,tensup(3,:),'--k');
hold on;
plot(tup(end)+t,tens(1,:),'-b',tup(end)+t,tens(2,:),...
     '-r',tup(end)+t,tens(3,:),'-k');
plot(tup(end)+t(end)+tdown,tensdown(1,:),...
    '-.b',tup(end)+t(end)+tdown,tensdown(2,:),...          
    '-.r',tup(end)+t(end)+tdown,tensdown(3,:),'-.k');
grid on;
xlabel('XXX');
ylabel('YYY');
lh = legend('111111111','222222222','333333333','444444444',...
            '555555555','666666666','777777777','888888888',...
            '999999999');
lh.Location = 'northeastoutside';

saveas(gca,'tension_exemple_appli.svg');