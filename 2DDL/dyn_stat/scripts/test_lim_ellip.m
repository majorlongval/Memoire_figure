%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script aims to test a function that determines the possible values
% of amplitude given the geometric parameters, the  ellipse parameters
% and the wrench.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;
addpath('functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L      = 5;            % m distance between independant and dep. pulleys
l      = 0.8;          % m width of end-effector
cx     = 0;            % m X position of centre of mass
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
a      = 0;            % m major axis
b      = 20;           % m minor axis
theta  = 0;            % rad angle of ellips
phi    = 0;            % rad starting angle
g      = 9.81;         % m/s^2 gravitational acceleration
omegaf = sqrt(g/12);   % rad/s angular frequency of traj
T      = 14*pi/omegaf; % s period of transition
tors   = mC*g*[0, 1 , 0];   % [N, N, Nm] The wrench

ry = 8;
k = 1/16;
rx = ry*k;
phix = 0;
phiy = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trajectory calculation  %%%%%%%%%%%%%%%%%%%%%

% [pos, acc, t,phix,phiy,rx,ry] = ...
%     traj2Dabtheta([x0;y0],a,b,theta,phi,omegaf,T);
% k = rx/ry;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating the conditions  %%%%%%%%%%%%%%%%%

[omega, omegan, rlimx_1, rlimx_2, rlimx_3,...
    rlimy_1, rlimy_2, rlimy_3] = ...
    ellip_cond(L, l, cx, cy, ay, x0, y0,phix,phiy, k, tors,m,res);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm to get limit vector  %%%%%%%%%%%%%%

rlim_y = [rlimy_1; rlimy_2; rlimy_3];
min = 1;
for i =1:length(rlim_y(:,1))
    if rlim_y(i,1)<= rlim_y(min,1)
        min = i;
        min_rlim_y = rlim_y(i,1);
    end
end

rlim_y_fin = rlim_y(min,1);
verif1 = 1; verif2 = 1; verif3 = 1;
for i =2:length(omega)
    [verif1, verif2, verif3] = ...
        ellip_cond_verif(L, l, cx, cy, ay, x0, y0,...
        phix,phiy, k, tors,m, omega(i),rlim_y(min,i));
    if verif1 == 1 && verif2 == 1 && verif3 == 1
        rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
    elseif verif1 == 0
        min = 1;
        rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
    elseif verif2 == 0
        min = 2;
        rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
    elseif verif3 == 0
        min = 3;
        rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
    end
end

rlim_x_fin = rlim_y_fin*k;















%% Printing the graph of the curves in rw
fig3 = figure;
set(gcf,'color','w');
plot(rlimx_1/x0,omega/omegan,':k','LineWidth',2); hold on; grid on;
plot(rlimx_2/x0,omega/omegan,'--k','LineWidth',2);
plot(rlimx_3/x0,omega/omegan,'-.k','LineWidth',2);
fh1 = fill([0,rlim_y_fin*k,0]/x0,...
     [0,omega,omega(end)]/omegan,[0.8,0.8,0.8]...
     ); hold on;
fh1.FaceAlpha = 0.6;
fh1.EdgeAlpha = 0;
% %plot(rx/x0,omegaf/omegan,'*r');
axis([0 0.5 0 2]);
xlabel('S01','FontSize',16); ylabel('S02','FontSize',16);
legend('S03','S04','S05');
% title('S06');
% 
plot([rx rx]/x0,[0 2],'r');


fig4 = figure;
set(gcf,'color','w');
 hold on;
% plot(ry/L,omegaf/omegan,'*r');
plot(rlimy_1/L,omega/omegan,':k','LineWidth',2);  grid on; hold on;
plot(rlimy_2/L,omega/omegan,'--k','LineWidth',2);
plot(rlimy_3/L,omega/omegan,'-.k','LineWidth',2);
fh2 = fill([0,rlim_y_fin,0]/L,...
     [0,omega,omega(end)]/omegan,[0.8,0.8,0.8]...
     );
fh2.FaceAlpha = 0.6;
fh2.EdgeAlpha = 0;
axis([0 7 0 4]);
% hold off;
xlabel('S01','FontSize',16); ylabel('S02','FontSize',16);
legend('S03','S04','S05');
% title('S06');

plot([ry ry]/L,[0 4],'r');

saveas(gca,'croisement_exemple.svg');



%% Calculating the limits of the interval
% for i =1:3
%     Psi(i) = 
%     X = m^2*(k^2*Psi(i)
%    omega_lim(i,1) =  
% end




% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Algorithm to get limit vector  %%%%%%%%%%%%%%
% 
% rlim_y = [rlimy_1; rlimy_2; rlimy_3];
% min = 1;
% for i =1:length(rlim_y(:,1))
%     if rlim_y(i,1)<= rlim_y(min,1)
%         min = i;
%         min_rlim_y = rlim_y(i,1);
%     end
% end
% 
% rlim_y_fin = rlim_y(min,1);
% verif1 = 1; verif2 = 1; verif3 = 1;
% for i =2:length(omega)
%     [verif1, verif2, verif3] = ...
%         ellip_cond_verif(L, l, cx, cy, ay, x0, y0,...
%         phix,phiy, k, tors,m, omega(i),rlim_y(min,i));
%     if verif1 == 1 && verif2 == 1 && verif3 == 1
%         rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
%     elseif verif1 == 0
%         min = 1;
%         rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
%     elseif verif2 == 0
%         min = 2;
%         rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
%     elseif verif3 == 0
%         min = 3;
%         rlim_y_fin = [rlim_y_fin,rlim_y(min,i)];
%     end
% end
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% plotting stuff  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% geo = pack_geo(ay,cy,cx,L,l);
% 
% %%%%%%%%% Calculating the tension throught the trajectory
% for i =1:length(pos)
% tens(:,i) = tension(m,geo,pos(:,i),acc(:,i),tors');
% end
% 
% %%%%%%%%% Plotting the tension graph
% fig1 = figure;
% set(gcf,'color','w');
% plot(t(1:1000),tens(1,1:1000),':b','LineWidth',2);
% hold on;
% plot(t(1:1000),tens(2,1:1000),'--b','LineWidth',2);
% plot(t(1:1000),tens(3,1:1000),'-.b','LineWidth',2);
% plot(t(1001:2000),tens(1,1001:2000),':r','LineWidth',2);
% plot(t(2001:end),tens(1,2001:end),':g','LineWidth',2);
% plot(t(1:400:end),tens(1,1:400:end),'ok');
% label = {'1','2','3','4','5','6','7','8'};
% labelpoints(t(1:400:end),tens(1,1:400:end),label,'position','NE');
% plot(t(1:1000),tens(2,1:1000),'--b','LineWidth',2);
% plot(t(1001:2000),tens(2,1001:2000),'--r','LineWidth',2);
% plot(t(2001:end),tens(2,2001:end),'--g','LineWidth',2);
% plot(t(1:400:end),tens(2,1:400:end),'ok');
% labelpoints(t(1:400:end),tens(2,1:400:end),label,'position','NE');
% plot(t(1:1000),tens(3,1:1000),'-.b','LineWidth',2);
% plot(t(1001:2000),tens(3,1001:2000),'-.r','LineWidth',2);
% plot(t(2001:end),tens(3,2001:end),'-.g','LineWidth',2);
% plot(t(1:400:end),tens(3,1:400:end),'ok');
% labelpoints(t(1:400:end),tens(3,1:400:end),label,'position','NE');
% axis([0 t(end) 0 max(tens(:))]);
% xlabel('time [s]','FontSize',16);
% ylabel('Tension [N]','FontSize',16);
% legend('Cable 1','Cable 2','Cable 3','Location','northoutside');
% 
% grid on;
% 
% 
% %%%%%%%% PLotting the trajectory in xy space
% fig2 = figure;
% set(gcf,'color','w');
% p1 = plot(pos(1,1:1000),pos(2,1:1000),'-b','LineWidth',2);
% hold on;
% p2 = plot(pos(1,1001:2000),pos(2,1001:2000),'-r','LineWidth',2);
% p3 = plot(pos(1,2001:end),pos(2,2001:end),'-g','LineWidth',2);
% p4 = plot([x0-rx x0+rx],[y0-ry y0-ry],'-m','LineWidth',2);
% plot([x0-rx x0+rx],[y0+ry y0+ry],'-m','LineWidth',2);
% plot([x0-rx x0-rx],[y0-ry y0+ry],'-m','LineWidth',2);
% plot([x0+rx x0+rx],[y0-ry y0+ry],'-m','LineWidth',2);
% plot(pos(1,1:400:end),pos(2,1:400:end),'ok');
% labelpoints(pos(1,1:400:end),pos(2,1:400:end),label,'position','NE');
% 
% xlabel('x [m]','FontSize',16);
% ylabel('y [m]','FontSize',16);
% legend([p1,p2,p3,p4],...
%     {'Rest to dynamic','Dynamic','Dynamic to rest','Bounding box'}...
%     ,'Location','northoutside');
% camroll(-90);
% axis([x0-rx-3 x0+rx+3 y0-ry-3 y0+ry+3]);
% ax = gca;
% ax.YAxisLocation ='right';
% grid on;
% 
% %%%%%%%%% Tracing the static workspace
% xSW = linspace(0,30,1000);
% [ymin,ymax] = y_min_stat_tors(geo,xSW,tors,m);
% p5 = plot(xSW,ymin,'-k', xSW,ymax,'-k',[0,0],[ymin(1),ymax(1)],'-k');
% 
% 
% 
% 
% 
% 
% 
% %% Printing the graph of the curves in rw
% fig3 = figure;
% set(gcf,'color','w');
% fill([0,rlim_y_fin*k,0]/x0,...
%      [0,omega,omega(end)]/omegan,[0.8,0.8,0.8]...
%      ); hold on;
% plot(rlimx_1/x0,omega/omegan,':k','LineWidth',2); hold on; grid on;
% plot(rlimx_2/x0,omega/omegan,'--k','LineWidth',2);
% plot(rlimx_3/x0,omega/omegan,'-.k','LineWidth',2);
% plot(rx/x0,omegaf/omegan,'*r');
% axis([0 2 0 2]);
% hold off;
% xlabel('S01','FontSize',16); ylabel('S02','FontSize',16);
% legend('S03','S04','S05');
% title('S06');
% 
% fig4 = figure;
% set(gcf,'color','w');
% fill([0,rlim_y_fin,0]/L,...
%      [0,omega,omega(end)]/omegan,[0.8,0.8,0.8]...
%      ); hold on;
% plot(ry/L,omegaf/omegan,'*r');
% plot(rlimy_1/L,omega/omegan,':k','LineWidth',2);  grid on;
% plot(rlimy_2/L,omega/omegan,'--k','LineWidth',2);
% plot(rlimy_3/L,omega/omegan,'-.k','LineWidth',2);
% axis([0 7 0 2]);
% hold off;
% xlabel('S01','FontSize',16); ylabel('S02','FontSize',16);
% legend('S03','S04','S05');
% title('S06');