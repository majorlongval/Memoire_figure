%% Calculating the boundaries of the wrench feasible workspace
% Author : Jordan M. Longval

%% Initialising
clear; clc; close all;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

% Setting the arbitrary values
R     = 1;         % m
r     = 0.2;       % m
alpha = 0;         % rad
rc    = 0.04;         % m
pc  = pi/6;      % rad
hc    = 0.5;         % m
% Setting the minimum and maximum heigth
zmin = 1;          % m
zmax = 20;         % m
% Range of values for the external wrench and the mass
Mx = [-2 2];     % Nm
My = [2 6];      % Nm
Mz = [1 1.5];    % Nm
fx = [-6 2];     % N
fy = [4 8];      % N
fz = 0;          % N
m  = 30;         % kg
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

zr = [zmin,zmax];


%% Calculating the WFW
[corners,volume,K] = WFW_Calculator(R,r,alpha,rc,pc,hc,m,wrench_comb,zr);


%% Plotting the whole WFW
fig1 = figure;
% plot3(points(:,1),points(:,2),points(:,3),'*r');
points = corners(:,4:6);
minx = 1.1*min(points(:,1));
maxx = 1.1*max(points(:,1));
miny = 1.1*min(points(:,2));
maxy = 1.1*max(points(:,2));
axis([minx maxx miny maxy 0.9*zmin 1.1*zmax])
hold on;
trh = trisurf(K,points(:,1),points(:,2),points(:,3));
trh.FaceColor = [0.2 0.2 1];
trh.FaceAlpha = 0.5;
trh.LineWidth = 2;
% Print projections for better view
proj_z = [points(:,1),points(:,2),ones(length(points(:,1)),1)*1.1*zmax];
Kpz = convhull(proj_z(:,1),proj_z(:,2));
plot3(points(Kpz,1),points(Kpz,2),1.1*zmax*ones(length(Kpz),1),'k--');
proj_x = [ones(length(points(:,1)),1)*maxx,points(:,2),points(:,3)];
Kpx = convhull(proj_x(:,2),proj_x(:,3));
plot3(maxx*ones(length(Kpx),1),points(Kpx,2),points(Kpx,3),'k--');
proj_y = [points(:,1),ones(length(points(:,1)),1)*miny,points(:,3)];
Kpy = convhull(proj_y(:,1),proj_y(:,3));
plot3(points(Kpy,1),miny*ones(length(Kpy),1),points(Kpy,3),'k--');

set (gca,'Ydir','reverse');
set (gca,'Zdir','reverse');
grid on;
xlabel('xlabel');
ylabel('ylabel');
zlabel('zlabel');
view(-43,15);
set(gcf,'renderer','Painters')
print -dsvg complet1.svg
%

%% Plotting the first step
planes = [];
for k =1:length(wrench_comb(:,1))
    planes = [planes;Calc_planes(alpha,R,r,rc,pc,hc,m,...
        wrench_comb(k,:))];
end
x = [minx maxx];
fig2 = figure;
for i =1:length(planes(:,1))
    y = (-1/planes(i,2))*(planes(i,1)*x+planes(i,3)*zmin+planes(i,4));
    lh(i)= plot3(x, y, [zmin,zmin],'-k');
    hold on;
end
for i =1:length(points(:,1))
    if points(i,3) == zmin
        plot3(points(i,1),points(i,2),points(i,3),'.r','MarkerSize',20);
    end
end

axis([minx maxx miny maxy 0.9*zmin 1.1*zmax])
set (gca,'Ydir','reverse');
set (gca,'Zdir','reverse');
grid on;
xlabel('xlabel');
ylabel('ylabel');
zlabel('zlabel');
view(-43,15);
set(gcf,'renderer','Painters')
print -dsvg top_lines.svg
%
%
%% Finding the next points
fig3 = figure;
proj_z = [points(:,1),points(:,2),ones(length(points(:,1)),1)*zmin];
Kpz = convhull(proj_z(:,1),proj_z(:,2));
plot3(points(Kpz,1),points(Kpz,2),zmin*ones(length(Kpz),1),'.r',...
    'MarkerSize',20);
ph = patch(points(Kpz,1),points(Kpz,2),zmin*ones(length(Kpz),1),[0.2,0.2,1]);
ph.FaceAlpha =0.5;
ph.LineWidth = 2;
view(3);
axis([minx maxx miny maxy 0.9*zmin 1.1*zmax])
set (gca,'Ydir','reverse');
set (gca,'Zdir','reverse');
grid on;
xlabel('xlabel');
ylabel('ylabel');
zlabel('zlabel');
view(-43,15);
set(gcf,'renderer','Painters')
hold on;
t = linspace(zmin-1,zmax+1,10);
pii = corners;
for i =1:length(pii(:,1))
    if pii(i,6) == zmin
        pxy_ini(:,:,i) = [planes(pii(i,1),1) planes(pii(i,1),2);...
            planes(pii(i,2),1) planes(pii(i,2),2)]\...
            -[planes(pii(i,1),3)*t+planes(pii(i,1),4);...
            planes(pii(i,2),3)*t+planes(pii(i,2),4)];
        lph(i) = plot3(pxy_ini(1,:,i),pxy_ini(2,:,i),t,'LineWidth',2);
    end
end

% Finding the points that are on each of the lines
set_of_lines = [];
for i =1:length(corners(:,1))
    if corners(i,6) == zmin
        set_of_lines = [set_of_lines;corners(i,1),corners(i,2),zmin];
    end
end

points_on_line = cell(length(set_of_lines(:,1)),1);
for i =1:length(set_of_lines(:,1))
    for j = 1:length(planes(:,1))
        temp = calc_plane_intersection(planes(set_of_lines(i,1),:),...
            planes(set_of_lines(i,2),:),...
            planes(j,:));
        if ~isempty(temp)
            temp2 = points_on_line{i};...
                temp2 = [temp2,[j;temp;set_of_lines(i,3)]];
            points_on_line{i} = temp2;
        end
    end
    p2p = points_on_line{i};
    p2p = p2p(:,p2p(4,:)>p2p(5,:));
    bph(i,j) = plot3(p2p(2,:),p2p(3,:),p2p(4,:),'.b','MarkerSize',20);
    jmin = 1;
    for j =1:length(p2p)
        if dist_2_points(corners(i,4:6),p2p(2:4,j))< dist_2_points(corners(i,4:6),p2p(2:4,jmin)) && p2p(4,j) > p2p(5,j)
            jmin = j;
        end
    end
    if p2p(4,jmin) > 0
    closest_point(:,i) = [corners(i,1);corners(i,2);p2p(:,jmin)];
    end
end
pause;
delete(bph);
delete(lph);
plot3(closest_point(4,:),closest_point(5,:),closest_point(6,:),'.r',...
     'MarkerSize',20);
%  
%% Repeting for the next step
set_of_lines = [];
for i =1:length(closest_point(1,:))
    set_of_lines = [set_of_lines;[closest_point(1,i),closest_point(3,i),closest_point(6,i)];...
                                 [closest_point(2,i),closest_point(3,i),closest_point(6,i)]];
end


points_on_line = cell(length(set_of_lines(:,1)),1);
for i =1:length(set_of_lines(:,1))
    for j = 1:length(planes(:,1))
        temp = calc_plane_intersection(planes(set_of_lines(i,1),:),...
            planes(set_of_lines(i,2),:),...
            planes(j,:));
        if ~isempty(temp)
            temp2 = points_on_line{i};...
                temp2 = [temp2,[j;temp;set_of_lines(i,3)]];
            points_on_line{i} = temp2;
        end
    end
    p2p = points_on_line{i};
    p2p = p2p(:,p2p(4,:)>p2p(5,:));
    bph(i,j) = plot3(p2p(2,:),p2p(3,:),p2p(4,:),'.b','MarkerSize',20);
    jmin = 1;
    for j =1:length(p2p)
        if dist_2_points(corners(i,4:6),p2p(2:4,j))< dist_2_points(corners(i,4:6),p2p(2:4,jmin)) && p2p(4,j) > p2p(5,j)
            jmin = j;
        end
    end
    if p2p(4,jmin) > 0
    closest_point(:,i) = [corners(i,1);corners(i,2);p2p(:,jmin)];
    end
end
pause;
delete(bph);
delete(lph);
plot3(closest_point(4,:),closest_point(5,:),closest_point(6,:),'.r',...
     'MarkerSize',20);

 
 
 

function dist = dist_2_points(p1,p2)
    dist = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2+(p1(3)-p2(3))^2);
end

function X = calc_plane_intersection(p1,p2,p3)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
M =[p1(1),p1(2),p1(3);...
    p2(1),p2(2),p2(3);...
    p3(1),p3(2),p3(3)];
if rank(M) == 3
    d = -[p1(4);p2(4);p3(4)];
    X = M\d;
else
    X = [];
    
    
    
end


end
