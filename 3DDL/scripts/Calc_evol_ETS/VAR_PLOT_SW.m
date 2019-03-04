%% VAR PLOT SW plots the evolution of the SW area over multiple values of rc(j)
clear all;clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
% Declaring the variables

R = 1;
r = 0.2;
alpha = (pi/2)-0.1;

hc = 0;
rc = linspace(0,2*r/3,11);
pc = 0;

sa = sin(alpha);
s2a = sin(2*alpha);
ca = cos(alpha);
c2a = cos(2*alpha);
spc = sin(pc);
cpc = cos(pc);

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
circle(0,0,R);
xlabel('x'); ylabel('y');
plot(R1(1),R1(2),'*r');
plot(R2(1),R2(2),'*r');
plot(R3(1),R3(2),'*r');
fillcolor = [1 0 0;1 0 0; 1 0 0];
    

for j =1:length(rc)

% Parameters of the linear functions
cx(1,1) = 2*r*ca^2+rc(j)*spc;
cx(1,2) = 2*r*ca^2-rc(j)*spc;
cx(2,1) = -r*ca^2+rc(j)*spc-sqrt(3)*r*s2a/2;
cx(2,2) = -r*ca^2-rc(j)*spc-sqrt(3)*r*s2a/2;
cx(3,1) = -r*ca^2+rc(j)*spc+sqrt(3)*r*s2a/2;
cx(3,2) = -r*ca^2-rc(j)*spc+sqrt(3)*r*s2a/2;

cy(1,1) = r*s2a-rc(j)*cpc;
cy(1,2) = r*s2a+rc(j)*cpc;
cy(2,1) = sqrt(3)*r*ca^2-r*s2a/2-rc(j)*cpc;
cy(2,2) = sqrt(3)*r*ca^2-r*s2a/2+rc(j)*cpc;
cy(3,1) = -sqrt(3)*r*ca^2-r*s2a/2-rc(j)*cpc;
cy(3,2) = -sqrt(3)*r*ca^2-r*s2a/2+rc(j)*cpc;

ck(1,1) = R*r*ca+2*R*rc(j)*ca*spc;
ck(1,2) = R*r*ca-2*R*rc(j)*ca*spc;
ck(2,1) = R*r*ca-R*rc(j)*ca*spc-sqrt(3)*R*rc(j)*ca*cpc;
ck(2,2) = R*r*ca+R*rc(j)*ca*spc+sqrt(3)*R*rc(j)*ca*cpc;
ck(3,1) = R*r*ca-R*rc(j)*ca*spc+sqrt(3)*R*rc(j)*ca*cpc;
ck(3,2) = R*r*ca+R*rc(j)*ca*spc-sqrt(3)*R*rc(j)*ca*cpc;

% setting the dimensions of the image



p.c(:,1) = [cx(1,1) cy(1,1);cx(1,2) cy(1,2)]\-[ck(1,1);ck(1,2)];
p.c(:,2) = [cx(1,1) cy(1,1);cx(2,1) cy(2,1)]\-[ck(1,1);ck(2,1)];
p.c(:,3) = [cx(1,1) cy(1,1);cx(2,2) cy(2,2)]\-[ck(1,1);ck(2,2)];
p.c(:,4) = [cx(1,1) cy(1,1);cx(3,1) cy(3,1)]\-[ck(1,1);ck(3,1)];
p.c(:,5) = [cx(1,1) cy(1,1);cx(3,2) cy(3,2)]\-[ck(1,1);ck(3,2)];
p.c(:,6) = [cx(1,2) cy(1,2);cx(2,1) cy(2,1)]\-[ck(1,2);ck(2,1)];
p.c(:,7) = [cx(1,2) cy(1,2);cx(2,2) cy(2,2)]\-[ck(1,2);ck(2,2)];
p.c(:,8) = [cx(1,2) cy(1,2);cx(3,1) cy(3,1)]\-[ck(1,2);ck(3,1)];
p.c(:,9) = [cx(1,2) cy(1,2);cx(3,2) cy(3,2)]\-[ck(1,2);ck(3,2)];
p.c(:,10) = [cx(2,1) cy(2,1);cx(2,2) cy(2,2)]\-[ck(2,1);ck(2,2)];
p.c(:,11) = [cx(2,1) cy(2,1);cx(3,1) cy(3,1)]\-[ck(2,1);ck(3,1)];
p.c(:,12) = [cx(2,1) cy(2,1);cx(3,2) cy(3,2)]\-[ck(2,1);ck(3,2)];
p.c(:,13) = [cx(2,2) cy(2,2);cx(3,1) cy(3,1)]\-[ck(2,2);ck(3,1)];
p.c(:,14) = [cx(2,2) cy(2,2);cx(3,2) cy(3,2)]\-[ck(2,2);ck(3,2)];
p.c(:,15) = [cx(3,1) cy(3,1);cx(3,2) cy(3,2)]\-[ck(3,1);ck(3,2)];


n =numel(fieldnames(p));
possib = nchoosek(1:15,2);
possib = [(1:length(possib(:,1)))',possib];

t = linspace(0,1,10);

lim_segs = [];
seg = struct;
for i =1:length(possib(:,1))
    [seg.mp(i,:),seg.len(i),seg.v(i,:)] = segment_calc(p.c(:,possib(i,2)),p.c(:,possib(i,3)));
    seg.p1(:,i) = p.c(:,possib(i,2));
    seg.p2(:,i) = p.c(:,possib(i,3));
    return_bool(i) = check_cond(seg.p1(:,i),seg.p2(:,i),R,r,rc(j),alpha,pc);
    if return_bool(i) == 1
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
fh = fill(seg(:,1),seg(:,2),fillcolor(mod(j,3)+1,:));
fh.EdgeColor = 'k';
fh.FaceAlpha = 1;
hold on;
end

saveas(gcf,'var_area_1','svg')






function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'-k');
end
