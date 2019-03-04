%% Script to draw the ETS as a function of the parameters
clear all;close all;clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
% Declaring the variables

R = 1;
r = 0.2;
alpha = pi/3;

hc = 0;
rc = 0.04;
pc = pi/3;

sa = sin(alpha);
s2a = sin(2*alpha);
ca = cos(alpha);
c2a = cos(2*alpha);
spc = sin(pc);
cpc = cos(pc);

% Parameters of the linear functions
cx(1,1) = 2*r*ca^2+rc*spc;
cx(1,2) = 2*r*ca^2-rc*spc;
cx(2,1) = -r*ca^2+rc*spc-sqrt(3)*r*s2a/2;
cx(2,2) = -r*ca^2-rc*spc-sqrt(3)*r*s2a/2;
cx(3,1) = -r*ca^2+rc*spc+sqrt(3)*r*s2a/2;
cx(3,2) = -r*ca^2-rc*spc+sqrt(3)*r*s2a/2;

cy(1,1) = r*s2a-rc*cpc;
cy(1,2) = r*s2a+rc*cpc;
cy(2,1) = sqrt(3)*r*ca^2-r*s2a/2-rc*cpc;
cy(2,2) = sqrt(3)*r*ca^2-r*s2a/2+rc*cpc;
cy(3,1) = -sqrt(3)*r*ca^2-r*s2a/2-rc*cpc;
cy(3,2) = -sqrt(3)*r*ca^2-r*s2a/2+rc*cpc;

ck(1,1) = R*r*ca+2*R*rc*ca*spc;
ck(1,2) = R*r*ca-2*R*rc*ca*spc;
ck(2,1) = R*r*ca-R*rc*ca*spc-sqrt(3)*R*rc*ca*cpc;
ck(2,2) = R*r*ca+R*rc*ca*spc+sqrt(3)*R*rc*ca*cpc;
ck(3,1) = R*r*ca-R*rc*ca*spc+sqrt(3)*R*rc*ca*cpc;
ck(3,2) = R*r*ca+R*rc*ca*spc-sqrt(3)*R*rc*ca*cpc;

% setting the dimensions of the image

xmax = 4.5;
xmin = -4.5;

x = linspace(xmin,xmax,1000);

Qs = [cos(2*pi/3) -sin(2*pi/3) 0;
    sin(2*pi/3)  cos(2*pi/3) 0;
    0                     0  1];
R1 = R*[cos(alpha);sin(alpha);0];
R2 = Qs*R1;
R3 = Qs*R2;
% 
figure;
hold on; grid on;
axis([-1.5 1.5 -1.5 1.5]);
axis square;
set(gca,'Ydir','reverse')
% plot(R1(1),R1(2),'*r');
% plot(R2(1),R2(2),'*r');
% plot(R3(1),R3(2),'*r');
circle(0,0,R);
xlabel('x'); ylabel('y');

for i =1:3
    for j =1:2
        if cy(i,j) ~= 0
            y(i,:,j) = (-cx(i,j)*x-ck(i,j))/cy(i,j);
            h_y(i,j) = plot(x,y(i,:,j));
        else
        h_y(i,j) = plot([-ck(i,j)/cx(i,j) -ck(i,j)/cx(i,j)],[-10 10]);
        end
    end
end

h_y(1,1).Color = 'r'; h_y(1,1).LineStyle = '-';
h_y(1,2).Color = 'r'; h_y(1,2).LineStyle = '--';
h_y(2,1).Color = 'b'; h_y(2,1).LineStyle = '-';
h_y(2,2).Color = 'b'; h_y(2,2).LineStyle = '--';
h_y(3,1).Color = 'g'; h_y(3,1).LineStyle = '-';
h_y(3,2).Color = 'g'; h_y(3,2).LineStyle = '--';

%delete(h_y(2,1));
% plot([R2(1) R3(1)],[R2(2) R3(2)],'-.r');
% plot([R1(1) R2(1)],[R1(2) R2(2)],'-.g');
% plot([R1(1) R3(1)],[R1(2) R3(2)],'-.b');

%Calcul des points d'intersection

p.c(:,1) = [cx(1,1) cy(1,1);cx(1,2) cy(1,2)]\-[ck(1,1);ck(1,2)];
p.c(:,2) = [cx(1,1) cy(1,1);cx(2,1) cy(2,1)]\-[ck(1,1);ck(2,1)];
p.c(:,3) = [cx(1,1) cy(1,1);cx(2,2) cy(2,2)]\-[ck(1,1);ck(2,2)];
p.c(:,4) = [cx(1,1) cy(1,1);cx(3,1) cy(3,1)]\-[ck(1,1);ck(3,1)];
p.c(:,5) = [cx(1,1) cy(1,1);cx(3,2) cy(3,2)]\-[ck(1,1);ck(3,2)];
p.c(:,6) = [cx(1,2) cy(1,2);cx(2,1) cy(2,1)]\-[ck(1,2);ck(2,1)];
%plot(pi.c(1,6),pi.c(2,6),'or');
p.c(:,7) = [cx(1,2) cy(1,2);cx(2,2) cy(2,2)]\-[ck(1,2);ck(2,2)];
p.c(:,8) = [cx(1,2) cy(1,2);cx(3,1) cy(3,1)]\-[ck(1,2);ck(3,1)];
%plot(pi.c(1,8),pi.c(2,8),'or');
p.c(:,9) = [cx(1,2) cy(1,2);cx(3,2) cy(3,2)]\-[ck(1,2);ck(3,2)];
p.c(:,10) = [cx(2,1) cy(2,1);cx(2,2) cy(2,2)]\-[ck(2,1);ck(2,2)];
p.c(:,11) = [cx(2,1) cy(2,1);cx(3,1) cy(3,1)]\-[ck(2,1);ck(3,1)];
p.c(:,12) = [cx(2,1) cy(2,1);cx(3,2) cy(3,2)]\-[ck(2,1);ck(3,2)];
%plot(pi.c(1,12),pi.c(2,12),'or');
p.c(:,13) = [cx(2,2) cy(2,2);cx(3,1) cy(3,1)]\-[ck(2,2);ck(3,1)];
p.c(:,14) = [cx(2,2) cy(2,2);cx(3,2) cy(3,2)]\-[ck(2,2);ck(3,2)];
p.c(:,15) = [cx(3,1) cy(3,1);cx(3,2) cy(3,2)]\-[ck(3,1);ck(3,2)];
%plot(pi.c(1,15),pi.c(2,15),'or');

%barycentre = (pi.c(:,6)+pi.c(:,8)+pi.c(:,12)+pi.c(:,15))/4;
%fill([pi.c(1,6),barycentre(1),pi.c(1,12)],[pi.c(2,6),barycentre(2),pi.c(2,12)],'r');
%fill([pi.c(1,12),barycentre(1),pi.c(1,15)],[pi.c(2,12),barycentre(2),pi.c(2,15)],'g');
%fill([pi.c(1,15),barycentre(1),pi.c(1,8)],[pi.c(2,15),barycentre(2),pi.c(2,8)],'m');
%fill([pi.c(1,8),barycentre(1),pi.c(1,6)],[pi.c(2,8),barycentre(2),pi.c(2,6)],'y');
%plot(barycentre(1),barycentre(2),'*b');
%plot([pi.c(1,6),barycentre(1)],[pi.c(2,6),barycentre(2)],'-b');
%plot([pi.c(1,12),barycentre(1)],[pi.c(2,12),barycentre(2)],'-b');
%plot([pi.c(1,8),barycentre(1)],[pi.c(2,8),barycentre(2)],'-b');
%plot([pi.c(1,15),barycentre(1)],[pi.c(2,15),barycentre(2)],'-b');
%plot([pi.c(1,6),pi.c(1,12)],[pi.c(2,6),pi.c(2,12)],'-b','LineWidth',2);
%plot([pi.c(1,6),pi.c(1,8)],[pi.c(2,6),pi.c(2,8)],'-b','LineWidth',2);
%plot([pi.c(1,8),pi.c(1,15)],[pi.c(2,8),pi.c(2,15)],'-b','LineWidth',2);
%plot([pi.c(1,15),pi.c(1,12)],[pi.c(2,15),pi.c(2,12)],'-b','LineWidth',2);
%lg_handle = legend([h_y(1,1) h_y(1,2) h_y(2,1) h_y(2,2) h_y(3,1) h_y(3,2)],...
%        {'S001','S002','S003','S004','S005','S006'});
% lg_handle.Location = 'eastoutside';

%saveas(gcf,'calc_area','svg')
n =numel(fieldnames(p));
possib = nchoosek(1:15,2);
% test = (1:length(possib(:,1)))';
possib = [(1:length(possib(:,1)))',possib];

t = linspace(0,1,10);
% Calcul des points milieu de chaque point d'intersection
lim_segs = [];
for i =1:length(possib(:,1))
    [seg.mp(i,:),seg.len(i),seg.v(i,:)] = segment_calc(p.c(:,possib(i,2)),p.c(:,possib(i,3)));
    seg.p1(:,i) = p.c(:,possib(i,2));
    seg.p2(:,i) = p.c(:,possib(i,3));
    return_bool(i) = check_cond(seg.p1(:,i),seg.p2(:,i),R,r,rc,alpha,pc);
    if return_bool(i) == 1
        %plot(seg.mp(i,1),seg.mp(i,2),'*b');
        lim_segs = [lim_segs,i];
    end
end

xseg = []; yseg = [];
for i =1:length(lim_segs)
    xseg =[xseg,p.c(1,possib(lim_segs(i),2)),p.c(1,possib(lim_segs(i),3))];
    yseg =[yseg,p.c(2,possib(lim_segs(i),2)),p.c(2,possib(lim_segs(i),3))];
end
seg1 = [xseg;yseg];
seg2 = unique(seg1','rows');
[seg_X,seg_Y] = vec_acw_order(seg2(:,1),seg2(:,2));
seg = [seg_X, seg_Y];
fh = fill(seg(:,1),seg(:,2),[0.6,0.6,1]);
fh.EdgeColor = 'none';
fh.FaceAlpha = 0.5;


 
area = polyarea(seg(:,1),seg(:,2));


fig2 = figure;
point_top = [seg,[0;0;0]];
point_fin = [seg,[10;10;10]];
fh1 = fill3([point_top(1,1) point_top(2,1)...
             point_fin(2,1) point_fin(1,1)],...
            [point_top(1,2) point_top(2,2)...
             point_fin(2,2) point_fin(1,2)],...
             [point_top(1,3) point_top(2,3)...
             point_fin(2,3) point_fin(1,3)],[0.6,0.6,1]);
hold on;
fh2 = fill3([point_top(1,1) point_top(3,1)...
             point_fin(3,1) point_fin(1,1)],...
            [point_top(1,2) point_top(3,2)...
             point_fin(3,2) point_fin(1,2)],...
             [point_top(1,3) point_top(3,3)...
             point_fin(3,3) point_fin(1,3)],[0.6,0.6,1]);
fh3 = fill3([point_top(2,1) point_top(3,1)...
             point_fin(3,1) point_fin(2,1)],...
            [point_top(2,2) point_top(3,2)...
             point_fin(3,2) point_fin(2,2)],...
             [point_top(2,3) point_top(3,3)...
             point_fin(3,3) point_fin(2,3)],[0.6,0.6,1]);
fh1.FaceAlpha = 0.5;
fh2.FaceAlpha = 0.5;
fh3.FaceAlpha = 0.5;

axis2 = gca;
set(axis2,'ZDir','reverse');
set(axis2,'YDir','reverse');
grid on;
circle3h = circle3(0,0,0,R);
plot3(R1(1),R1(2),R1(3),'*r');
plot3(R2(1),R2(2),R2(3),'*b');
plot3(R3(1),R3(2),R3(3),'*g');
xlabel('XLABEL');
ylabel('YLABEL');
zlabel('ZLABEL');
title({'TITLE1','TITLE2'});
%print('-painters','ETS_TOT','-dsvg');
function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'-k');
end

function h = circle3(x,y,z,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    zunit = z*ones(length(th));
    h = plot3(xunit, yunit,zunit,'-k');
end














