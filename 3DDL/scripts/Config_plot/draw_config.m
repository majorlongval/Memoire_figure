%% Script to draw config
clear all; close all; clc;

% Declaring the position of the end-effector
x = 0; y = 0; z = 2;
p = [x;y;z];
% Declaring the construction values
R = 1;
r = 0.2;
alpha = 1.0747;

% Declaring the vectors
a11 = r*[0; 1; 0];
a12 = -a11;
ct = cos(2*pi/3); st = sin(2*pi/3);
Qs = [ct -st 0;st ct 0;0 0 1];
a21 = Qs*a11;
a22 = -a21;
a31 = Qs*a21;
a32 = -a31;
R1 = R*[cos(alpha);sin(alpha);0];
R2 = Qs*R1;
R3 = Qs*R2;

t = linspace(0,1,100); % parameter for length of cable

% Defining the lines
line11 = R1+a11 + (p-R1).*t;
line12 = R1+a12 + (p-R1).*t;
line21 = R2+a21 + (p-R2).*t;
line22 = R2+a22 + (p-R2).*t;
line31 = R3+a31 + (p-R3).*t;
line32 = R3+a32 + (p-R3).*t;

figure
p = plot3(x,y,z,'*r');
hold on; grid on;
h_base = circle(0,0,0,R);
h_eff = circle(x,y,z,r);
% Drawing the lines
l11 = plot3(line11(1,:),line11(2,:),line11(3,:),'Color',...
      [0 0 1]);
l12 = plot3(line12(1,:),line12(2,:),line12(3,:),'Color',...
      [0 0 0.5]);
l21 = plot3(line21(1,:),line21(2,:),line21(3,:),'Color',...
      [1 0 0]);
l22 = plot3(line22(1,:),line22(2,:),line22(3,:),'Color',...
      [0.5 0 0]);
l31 = plot3(line31(1,:),line31(2,:),line31(3,:),'Color',...
      [0 1 0]);
l32 = plot3(line32(1,:),line32(2,:),line32(3,:),'Color',...
      [0 0.5 0]);

%camorbit(180,0,'data',[1 0 0]);
camup([0 0 -1]);
cameratoolbar('Show');
xlabel('x');
ylabel('y');
zlabel('z');
axis equal ;
leg = legend([p,l11,l12,l21,l22,l31,l32],{'P','l11long','l12long','l21long','l22long','l31long','l32long'});
leg.FontSize = 14;


% saveas(gcf,'alpha_p2_crois','svg');

function h = circle(x,y,z,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
zunit = z*ones(length(th));
h = plot3(xunit, yunit,zunit,'-k');
end

