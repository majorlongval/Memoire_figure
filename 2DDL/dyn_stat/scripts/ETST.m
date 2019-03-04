%% Script to determine the ETST

clear all; close all; clc;

addpath('functions');

%% Input values

% weight and wrench
g = 9.81;
mR = 2;
mC = 1;

% Geometric parameters
% The unpack and pack is really a bad idea
[ay, cy, cx, L, l,A, B, C] = unpack_geo(pack_geo(0.3,-0.1,0.2,5,0.8));

dy = -0.2;

Tors1 = [mC*g 0 dy*mC*g];

Tors2 = [0 0 0];


D(1) = 0;
D(2) = 2*L+ay;
D(3) = -D(2);

pente = Tors1(2)/(mR*g+Tors1(1));

b(1:3) = (C(1:3).*Tors1(2)-A(1:3).*(mR*g+Tors1(1))-D(1:3).*Tors1(3))./...
         (B(1:3).*(mR*g+Tors1(1)));

x = linspace(0, 5,1000);

y(:,1) = pente*x+b(1);
y(:,2) = pente*x+b(2);
y(:,3) = pente*x+b(3);

y2(:,1) = ones(1,1000)*-A(1)/B(1);
y2(:,2) = ones(1,1000)*-A(2)/B(2);
y2(:,3) = ones(1,1000)*-A(3)/B(3);

figure


% p0 = fill([x,x(end),fliplr(x),x(1)],[y(:,3)',...
%     y(end,1),fliplr(y(:,1)'),y(1,1)],[0.8,0.8,0.8],'EdgeColor','none');
hold on;
p1 = plot(x,y(:,1),'-b','LineWidth',2);
p2 = plot(x,y(:,2),'-r','LineWidth',2);
p3 = plot(x,y(:,3),'-k','LineWidth',2);




camroll(-90);

axis([0 5 -10 10]);
ax=gca;
ax.YAxisLocation = 'Right';
xlabel('S01');
ylabel('S02');
title('S03');

% l1 = legend([p0, p1,p2,p3,p4,p5,p6],{'S04','S05','S06','S07','S08',...
%                                      'S09','S10'});
% l1.Location =  'southwest';



Tors2 = [0 0 0];

hold on;
D(1) = 0;
D(2) = 2*L+ay;
D(3) = -D(2);

pente = Tors2(2)/(mR*g+Tors2(1));

b(1:3) = (C(1:3).*Tors2(2)-A(1:3).*(mR*g+Tors2(1))-D(1:3).*Tors2(3))./...
         (B(1:3).*(mR*g+Tors2(1)));

x = linspace(0, 5,1000);

y(:,1) = pente*x+b(1);
y(:,2) = pente*x+b(2);
y(:,3) = pente*x+b(3);

y2(:,1) = ones(1,1000)*-A(1)/B(1);
y2(:,2) = ones(1,1000)*-A(2)/B(2);
y2(:,3) = ones(1,1000)*-A(3)/B(3);



% p0 = fill([x,x(end),fliplr(x),x(1)],[y(:,3)',...
%     y(end,1),fliplr(y(:,1)'),y(1,1)],[0.8,0.8,0.8],'EdgeColor','none');
hold on;
p4 = plot(x,y(:,1),'--b','LineWidth',2);
p5 = plot(x,y(:,2),'--r','LineWidth',2);
p6 = plot(x,y(:,3),'--k','LineWidth',2);
l1 = legend([p1 p3 p4 p6],{'1','2','3','4'});
l1.Location =  'southwest';


plot(5,0,'*r');
plot(4,0,'*b');
plot(3,0.5,'*g');
plot(2.5,1,'*m');

grid on;


saveas(gca,'EXEMPLE_APPLI.svg');
%print('ETST','-depsc2','-painters');
