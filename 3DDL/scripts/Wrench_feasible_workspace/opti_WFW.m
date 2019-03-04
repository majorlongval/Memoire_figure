%% Script to omptimise geometry of EE for volume of WFW
clear; close all; clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
set(groot, 'defaultFigureRenderer', 'painters')
% Setting the arbitrary values
R     = 1;         % m
r     = 0.3;       % m
alpha = pi/6;         % rad
rc    = 0;         % m
phic  = pi/6;      % rad
hc    = 0.00001;         % m


% Setting the minimum and maximum heigth
zmin = 1;          % m
zmax = 10;         % m
zr = [zmin,zmax];

% Range of values for the external wrench and the mass
Mx = 2; % Nm
My = -1; % Nm
Mz = 2; % Nm
fx = -20; % N
fy = -10; % N
fz = 0;      % N
m  = 10;     % kg
% Defining all the possible combinations of wrenches
wrench_comb = [];
for a = 1:length(fx)
    for b =1:length(fy)
        for c = 1:length(fz)
            for d = 1:length(Mx)
                for e = 1:length(My)
                    for f = 1:length(Mz)
                        wrench_comb = [wrench_comb;...
                            fx(a),fy(b),fz(c),Mx(d),My(e),Mz(f)];
                    end
                end
            end
        end
    end
end

wrench_comb2 = [0.0001,0,0,0,0,0];
fun = @(x)-WFW_Calculator(R,r,x(1),x(2),x(3),x(4),m,wrench_comb,zr);
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
lb = [0,0,0,-r];
ub = [pi/2,r,pi/3,r];
x0 = (lb+ub)/2;
options = optimoptions('fmincon','Display','iter')
%options = optimoptions('fmincon','Algorithm','sqp');
[x,opti_vol] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
alpha_sol = x(1);
rc_sol = x(2);
phic_sol = x(3);
h_sol = x(4);
% 
opti_vol = opti_vol*-1;
% 
[vol,K,vertices] = WFW_Calculator(R,r,alpha,rc,phic,hc,m,wrench_comb,zr);
[vol2,K2,vertices2] = WFW_Calculator(R,r,alpha_sol,rc_sol,phic_sol,h_sol,m,wrench_comb,zr);
% % Printing the robot setup
[fig1_h,fig2_h] = Robot_SETUP(R,r,alpha,rc,phic,hc);
% Print projections for better view


xlabel(fig1_h,'XXX');
ylabel(fig1_h,'YYY');
zlabel(fig1_h,'ZZZ');
axis(fig1_h,[-1.1,1.1,-1.1,1.1]);
axis(fig1_h,'square');
set(fig1_h,'YDir','reverse');
set(fig1_h,'ZDir','reverse');

axis(fig2_h,'square');
axis(fig2_h,[-1.5,1.5,-1.5,1.5,0,2.5]);
set(fig2_h,'YDir','reverse');
set(fig2_h,'ZDir','reverse');
xlabel(fig2_h,'XXX');
ylabel(fig2_h,'YYY');
zlabel(fig2_h,'ZZZ');

figure
tsh = trisurf(K,vertices(:,4),vertices(:,5),vertices(:,6));
tsh.FaceColor = 'r';
tsh.FaceAlpha = 0.3;
axis('square');
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
axis([-2 2 -2 2 0 10]);
xlabel('XXX');
ylabel('YYY');
zlabel('ZZZ');

hold on;
points = vertices(:,4:6);
minx = 2*min(points(:,1));
maxx = 2*max(points(:,1));
miny = 2*min(points(:,2));
maxy = 2*max(points(:,2));
minz = 0.9*zmin;
maxz = 1.1*zmax;
axis([minx maxx miny maxy minz maxz]);
proj_z = [points(:,1),points(:,2),ones(length(points(:,1)),1)*1.1*zmax];
Kpz = convhull(proj_z(:,1),proj_z(:,2));
plot3(points(Kpz,1),points(Kpz,2),1.1*zmax*ones(length(Kpz),1),'k--');
proj_x = [ones(length(points(:,1)),1)*maxx,points(:,2),points(:,3)];
Kpx = convhull(proj_x(:,2),proj_x(:,3));
plot3(maxx*ones(length(Kpx),1),points(Kpx,2),points(Kpx,3),'k--');
proj_y = [points(:,1),ones(length(points(:,1)),1)*miny,points(:,3)];
Kpy = convhull(proj_y(:,1),proj_y(:,3));
plot3(points(Kpy,1),miny*ones(length(Kpy),1),points(Kpy,3),'k--');

saveas(gca,'fig3','epsc');
saveas(fig1_h,'fig1','epsc');
saveas(fig2_h,'fig2','epsc');






[fig3_h,fig4_h] = Robot_SETUP(R,r,alpha_sol,rc_sol,phic_sol,h_sol);
axis(fig3_h,[-1.1,1.1,-1.1,1.1]);
axis(fig3_h,'square');
set(fig3_h,'YDir','reverse');
set(fig3_h,'ZDir','reverse');
xlabel(fig3_h,'XXX');
ylabel(fig3_h,'YYY');
zlabel(fig3_h,'ZZZ');
axis(fig4_h,'square');
set(fig4_h,'YDir','reverse');
set(fig4_h,'ZDir','reverse');
axis('square');
axis(fig4_h,[-1.5,1.5,-1.5,1.5,0,2.5]);
xlabel(fig4_h,'XXX');
ylabel(fig4_h,'YYY');
zlabel(fig4_h,'ZZZ');

figure;
tsh2 = trisurf(K2,vertices2(:,4),vertices2(:,5),vertices2(:,6));
tsh2.FaceColor = 'b';
tsh2.FaceAlpha = 0.3;
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
axis('square');
xlabel('XXX');
ylabel('YYY');
zlabel('ZZZ');

points = vertices2(:,4:6);
minx = 2*min(points(:,1));
maxx = 2*max(points(:,1));
miny = 2*min(points(:,2));
maxy = 2*max(points(:,2));
minz = 0.9*zmin;
maxz = 1.1*zmax;
axis([minx maxx miny maxy minz maxz]);
hold on;
proj_z = [points(:,1),points(:,2),ones(length(points(:,1)),1)*1.1*zmax];
Kpz = convhull(proj_z(:,1),proj_z(:,2));
plot3(points(Kpz,1),points(Kpz,2),1.1*zmax*ones(length(Kpz),1),'k--');
proj_x = [ones(length(points(:,1)),1)*maxx,points(:,2),points(:,3)];
Kpx = convhull(proj_x(:,2),proj_x(:,3));
plot3(maxx*ones(length(Kpx),1),points(Kpx,2),points(Kpx,3),'k--');
proj_y = [points(:,1),ones(length(points(:,1)),1)*miny,points(:,3)];
Kpy = convhull(proj_y(:,1),proj_y(:,3));
plot3(points(Kpy,1),miny*ones(length(Kpy),1),points(Kpy,3),'k--');


saveas(gca,'test3','epsc');
% Checking the tension at the points
figure;
for i =1:length(vertices(:,1))
    tens(:,i) = calc_tension_traj3D(m,R,r,alpha,rc,phic,hc,wrench_comb,vertices(i,4:6)',[0;0;0]);
    plot(i*ones(length(tens(:,i))),tens(:,i),'*r');
    hold on;
end
grid on;

figure
for i =1:length(vertices(:,1))
    tens2(:,i) = calc_tension_traj3D(m,R,r,alpha_sol,rc_sol,phic_sol,h_sol,wrench_comb,vertices2(i,4:6)',[0;0;0]);
    plot(i*ones(length(tens2(:,i))),tens2(:,i),'*r');
    hold on;
end
grid on;

saveas(fig3_h,'fig4','epsc');
saveas(fig4_h,'fig5','epsc');