%% Script to plot the end-effector for a given configuration
clear all; close all;clc;
r =0.2;
r_cm = 0.01;
rc = 0.1;
pc = pi/8;
p = [rc*cos(pc),rc*sin(pc)];
h_cm = plot_cm_logo(p,r_cm);

h = circle(0,0,r)

a11 = [0;r];
a12 = -a11;
Qs = [cos(2*pi/3) -sin(2*pi/3);
      sin(2*pi/3)  cos(2*pi/3)];
a21 = Qs*a11;
a22 = -a21;
a31 = Qs*a21;
a32 = -a31;
plot(0,0,'xk');
plot(a11(1),a11(2),'xk');
plot(a12(1),a12(2),'xk');
plot(a21(1),a21(2),'xk');
plot(a22(1),a22(2),'xk');
plot(a31(1),a31(2),'xk');
plot(a32(1),a32(2),'xk');
plot([0 r],[0 0],'--k');
plot([0 r*cos(pc)],[0 r*sin(pc)],'--k');
plot([0 rc*cos(pc)],[0 rc*sin(pc)],'-k');
plot(((r+rc)/2)*cos(0:pi/200:pc),((r+rc)/2)*sin(0:pi/200:pc),'-k');
Qs = [cos(2*pi/3) -sin(2*pi/3);
      sin(2*pi/3)  cos(2*pi/3)];
      
axis('equal');
axis('square');
grid on;
set(gca,'Ydir','reverse')

















function h_cm = plot_cm_logo(p,r_cm)
ang = 0:2*pi/500:2*pi;
c = [0 0 0;
     1 1 1;
     0 0 0;
     1 1 1];
for i =1:4
    deb = (i-1)*length(ang)/4;
    fin = i*length(ang)/4;
    cmx = p(1);
    cmx = [cmx,p(1)+r_cm.*cos(ang(deb:fin))];
    cmx = [cmx,p(1)];
    cmy = p(2);
    cmy = [cmy,p(2)+r_cm.*sin(ang(deb:fin))];
    cmy = [cmy,p(2)];
    
    h_cm(i) = patch(cmx,cmy,c(i,:));
    hold on;
end
end


function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'-k');
end