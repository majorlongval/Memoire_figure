% TEST_TRAJ_TENS
clear all; close all; clc;

addpath('functions');

% GEOMETRY
m =1;
g = 9.81;
L = 5;
l = 1;
ay = 0.2;
cy = -0.2;
cx = 0.1;
x0 = 8;
y0 = 2;
% TRAJECTORY
a = 5;
b = 10;
theta = pi/12;
phi = 0;
omegaf = 0.9*sqrt(g/8);
T = 3.2*pi/omegaf;



[pos, acc, t,phix,phiy,rx,ry] = ...
    traj2Dabtheta([x0;y0],a,b,theta,phi,omegaf,T);

geo = pack_geo(ay,cy,cx,L,l);



% % Creation d'un torseur ressort
% kr = 5; %N/m
% pos_b = [10;0];
% l0 = norm([x0;y0]-pos_b);
% for i =1:length(pos(1,:))
%    er(:,i) =  pos(:,i)-pos_b;
%    lr = norm(er(:,i));
%    er(:,i) = er(:,i)/lr;
%    Tors(:,i) = -kr*(lr-l0)*[er(:,i);-(cx*er(2,i)-cy*er(1,i));];
% end
% 

Tors = [0;0;0];

for i =1:length(pos)
tens(:,i) = tension(m,geo,pos(:,i),acc(:,i),Tors);
end

fig1 = figure;
set(gcf,'color','w');
plot(t(1:1000),tens(1,1:1000),':b','LineWidth',2);
hold on;
plot(t(1:1000),tens(2,1:1000),'--b','LineWidth',2);
plot(t(1:1000),tens(3,1:1000),'-.b','LineWidth',2);
plot(t(1001:2000),tens(1,1001:2000),':r','LineWidth',2);
plot(t(2001:end),tens(1,2001:end),':g','LineWidth',2);
plot(t(1:400:end),tens(1,1:400:end),'ok');
label = {'1','2','3','4','5','6','7','8'};
labelpoints(t(1:400:end),tens(1,1:400:end),label,'position','NE');
plot(t(1:1000),tens(2,1:1000),'--b','LineWidth',2);
plot(t(1001:2000),tens(2,1001:2000),'--r','LineWidth',2);
plot(t(2001:end),tens(2,2001:end),'--g','LineWidth',2);
plot(t(1:400:end),tens(2,1:400:end),'ok');
labelpoints(t(1:400:end),tens(2,1:400:end),label,'position','NE');
plot(t(1:1000),tens(3,1:1000),'-.b','LineWidth',2);
plot(t(1001:2000),tens(3,1001:2000),'-.r','LineWidth',2);
plot(t(2001:end),tens(3,2001:end),'-.g','LineWidth',2);
plot(t(1:400:end),tens(3,1:400:end),'ok');
labelpoints(t(1:400:end),tens(3,1:400:end),label,'position','NE');
axis([0 t(end) 0 max(tens(:))]);
xlabel('xlabel','FontSize',16);
ylabel('ylabel','FontSize',16);
legend('Cable1','Cable2','Cable3','Location','northoutside');

grid on;

saveas(fig1,'tension.eps','epsc');

fig2 = figure;
set(gcf,'color','w');
p1 = plot(pos(1,1:1000),pos(2,1:1000),'--b','LineWidth',2);
hold on;
p2 = plot(pos(1,1001:2000),pos(2,1001:2000),'-r','LineWidth',2);
p3 = plot(pos(1,2001:end),pos(2,2001:end),'-.g','LineWidth',2);
% p4 = plot([x0-rx x0+rx],[y0-ry y0-ry],'-m','LineWidth',2);
% plot([x0-rx x0+rx],[y0+ry y0+ry],'-m','LineWidth',2);
% plot([x0-rx x0-rx],[y0-ry y0+ry],'-m','LineWidth',2);
% plot([x0+rx x0+rx],[y0-ry y0+ry],'-m','LineWidth',2);
p5 = plot([0 20],[y_min_stat(geo) y_min_stat(geo)],':k','LineWidth',2);
plot([0 20],[L L],':k','LineWidth',2);
plot([0 0],[y_min_stat(geo) L],':k','LineWidth',2);
plot(pos(1,1:400:end),pos(2,1:400:end),'ok');
labelpoints(pos(1,1:400:end),pos(2,1:400:end),label,'position','NE');

xlabel('xlabel','FontSize',16);
ylabel('ylabel','FontSize',16);

camroll(-90);
axis([x0-rx-6 x0+rx+3 y0-ry-3 y0+ry+3]);
ax = gca;
ax.YAxisLocation ='right';
grid on;
gridLegend([p1,p2,p3,p5],2,...
           {'legend1','legend2','legend3','legend5'},...
           'location','north','Fontsize',16);
       
saveas(fig2,'path_traj.eps','epsc');





%% Printing the graph of the curves in rw


k = rx/ry;

[omega, omegan, rlimx_1, rlimx_2, rlimx_3,...
    rlimy_1, rlimy_2, rlimy_3] = ...
    ellip_cond(L, l, cx, cy, ay, x0, y0, phix,phiy, k, Tors, m, 10000);

fig3 = figure;
set(gcf,'color','w');
plot(rlimx_1/x0,omega/omegan,':k','LineWidth',2); hold on; grid on;
plot(rlimx_2/x0,omega/omegan,'--k','LineWidth',2);
plot(rlimx_3/x0,omega/omegan,'-.k','LineWidth',2);
plot(rx/x0,omegaf/omegan,'*r');
axis([0 2 0 2]);
xticks([0 0.5 1 1.5 2]); yticks([0 0.5 1 1.5 2]);
hold off;
xlabel('S01','FontSize',16); ylabel('S02','FontSize',16);
l3 = legend('S0003','S0004','S0005');
l3.FontSize = 16;
title('S06');

saveas(fig3,'rxomega.eps','epsc');


% Finding a value for a minimum value of T
