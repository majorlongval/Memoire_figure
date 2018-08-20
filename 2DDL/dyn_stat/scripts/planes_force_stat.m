%% Planes of static forces
clear all;close all;clc
addpath('functions');
% input values
x  = 5;
cy = linspace(-1,1,50);
cx = 0;
ay = cy;
l  = 1;
L  = 4;
g = 9.81;

for j =1:length(cy)
    for i =1:length(ay)
        geo = pack_geo(ay(i),cy(j),cx,L,l);
        ymin = y_min_stat(geo);
        y = linspace(ymin,L,50);
        for k =1:length(y)
            [V(k),Me(k,:),fy(k,:)] = Effort_poly(x, y(k), geo);
        end
        V(find(isnan(V)))=[];
        V_m(i,j) = mean(V);
        Memax_m(i,j) = mean(Me(:,1));
        Memin_m(i,j) = mean(Me(:,2));
        Me_m(i,j) =    mean(Me(:,2)-Me(:,1));
        fy_m(i,j) =    mean(fy(:,2)-fy(:,1)); 
        fymax_m(i,j) = mean(fy(:,1));
        fymin_m(i,j) = mean(fy(:,2));
    end
    display(strcat('j=',num2str(j)));
end

figure
surf(ay,cy,V_m)

figure
surf(ay,cy,Me_m)

figure
surf(ay,cy,fy_m)
















% Derived values
% A1 = -L;
% A2 = ay*(cy-L)+2*L*cy-l*(ay+L);
% A3 = ay*(L-cy)-2*L*cy-l*(ay+L);
%
% B1 = 1;
% B2 = ay-l;
% B3 = -(ay+l);
%
% D1 = 0;
% D2 = 2*L+ay;
% D3 = -D2;
%
% C1 = 0;
% C2 = cx*D2;
% C3 = -C2;
%
%
% % plotting the planes
% [fy Me] = meshgrid(-10:0.5:10);
% %p1: m1*fx+n1*fy+p1*Me+q1=0
% m1 = A1+B1*y;
% n1 = C1+B1*x;
% p1 = D1;
% q1 = m1*g;
% fx1 = -(1/m1)*(n1*fy+p1*Me+q1);
%
% surf(fy,Me,fx1);
%
% %p2: m2*fx+n2*fy+p2*Me+q2=0
% m2 = A2+B2*y;
% n2 = C2+B2*y;
% p2 = D2;
% q2 = m2*g;
% fx2 = -(1/m2)*(n2*fy+p2*Me+q2);
% hold on;
% surf(fy,Me,fx2);
%
% %p2: m2*fx+n2*fy+p2*Me+q2=0
% m3 = A3+B3*y;
% n3 = C3+B3*x;
% p3 = D3;
% q3 = m3*g;
% fx3 = -(1/m3)*(n3*fy+p3*Me+q3);
% hold on;
% surf(fy,Me,fx3);
% xlabel('fy'); ylabel('Me'); zlabel('fx');
%
% % Trouver les valeurs fy et Me @ fx = 0
% Mat12 = [n1 p1; n2 p2];
% b12 = -g*[m1;m2];
% p12 = Mat12\b12;
% plot(p12(1),p12(2),'*r');
% u12 = [p12;0]-[0; 0; -g];
% plot3([0 p12(1)],[0 p12(2)],[-g 0],'-r','LineWidth',2);
%
% Mat13 = [n1 p1; n3 p3];
% b13 = -g*[m1;m3];
% p13 = Mat13\b13;
% plot(p13(1),p13(2),'*r');
% u13 = [p13;0]-[0; 0; -g];
% plot3([0 p13(1)],[0 p13(2)],[-g 0],'-r','LineWidth',2);
%
% Mat23 = [n2 p2; n3 p3];
% b23 = -g*[m2;m3];
% p23 = Mat23\b23;
% u23 = [p23;0]-[0; 0; -g];
% plot(p23(1),p23(2),'*r');
% plot3([0 p23(1)],[0 p23(2)],[-g 0],'-r','LineWidth',2);
%
% % Plotting extra lines to show tetrahedron
% plot3([p12(1) p23(1)],[p12(2) p23(2)],[0 0],'-r','LineWidth',2);
% plot3([p23(1) p13(1)],[p23(2) p13(2)],[0 0],'-r','LineWidth',2);
% plot3([p13(1) p12(1)],[p13(2) p12(2)],[0 0],'-r','LineWidth',2);
%
%
% % Calculating the values of range of Me at fy = fx = 0 and plot
% Me2 = -m2*g/p2;
% Me3 = -m3*g/p3;
% Mer = sort([-(A2+B2*y)*g/D2;-(A3+B3*y)*g/D3]);
% plot3([0 0],[Me2, Me3],[0 0],'-b','LineWidth',2);
%
%
% % Calculating the values of range of fy at Me = fx = 0 and plot
% fy1 = -m1*g/n1;
% fy2 = -m2*g/n2;
% fy3 = -m3*g/n3;
% fyr = [min_abs(fy2,fy3);fy1];
% plot3([fy1 fyr(1)],[0,0],[0 0],'-b','LineWidth',2);
%
% % Calculating the Volume
% V = abs(det([u12,u13,u23]));