%% Script pour tester l'équivalence entre les équations normale et simple
clear all; close all; clc;

%% Définition des paramètres normaux

% major and minor axis
a = 5; %major (m)
b = 8; %minor (m)

% frequency and starting angle
omega = pi;
phi   = -3*pi/2;

% rotation angles
alpha = 20; 
beta  = 5;
gamma = -210;

% rotation matrix of the active form in sequance Z1Y2X3
s1 = sin(alpha);
c1 = cos(alpha);
s2 = sin(beta);
c2 = cos(beta);
s3 = sin(gamma);
c3 = cos(gamma);
Q  = [c1*c2 c1*s2*s3-c3*s1 s1*s3+c1*c3*s2;...
      c2*s1 c1*c3+s1*s2*s3 c3*s1*s2-c1*s3;...
      -s2   c2*s3          c2*c3];

% centre of the ellipses;
p0 = [0; 0; 0];

%% tracer le graph normal;

% discrétisation du temps
t = linspace(0,2,100);

p(1:3,:) = p0 + Q*[a*cos(omega*t+phi);b*sin(omega*t+phi);zeros(1,length(t))];

comet3(p(1,:),p(2,:),p(3,:));
hold on;
plot3(p(1,1),p(2,1),p(3,1),'*r');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');


%% Définir les paramètres de la forme simplifiée
q11 = Q(1,1); q12 = Q(1,2); q13 = Q(1,3);
q21 = Q(2,1); q22 = Q(2,2); q23 = Q(2,3);
q31 = Q(3,1); q32 = Q(3,2); q33 = Q(3,3);

cp = cos(phi); sp = sin(phi);
MM = [q11*a*cp+q12*b*sp q12*b*cp-q11*a*sp;...
      q21*a*cp+q22*b*sp q22*b*cp-q21*a*sp;...
      q31*a*cp+q32*b*sp q32*b*cp-q31*a*sp];

px = atan2(MM(1,1),MM(1,2));
py = atan2(MM(2,1),MM(2,2));
pz = atan2(MM(3,1),MM(3,2));

rx =  ...%((MM(1,1)*sqrt(MM(1,1)^2))+(MM(1,2)*sqrt((MM(1,2)^2))))/...
      sqrt(MM(1,1)^2+MM(1,2)^2);
ry =  ...%((MM(2,1)*sqrt(MM(2,1)^2))+(MM(2,2)*sqrt(MM(2,2)^2)))/...
      sqrt(MM(2,1)^2+MM(2,2)^2);
rz =  ...%((MM(3,1)*sqrt(MM(3,1)^2))+(MM(3,2)*sqrt(MM(3,2)^2)))/...
      sqrt(MM(3,1)^2+MM(3,2)^2);
  
p2 = p0 + [rx*sin(omega*t+px);...
           ry*sin(omega*t+py);...
           rz*sin(omega*t+pz)];
       
comet3(p2(1,:),p2(2,:),p2(3,:));
plot3(p2(1,1),p2(2,1),p2(3,1),'*g');




%% Plotting a box rx,ry,rz
c1 = p0 + [rx; ry; rz];
c2 = p0 + [rx; ry; -rz];
c3 = p0 + [rx; -ry; rz];
c4 = p0 + [-rx; ry; rz];
c5 = p0 + [rx; -ry; -rz];
c6 = p0 + [-rx; ry; -rz];
c7 = p0 + [-rx; -ry; rz];
c8 = p0 + [-rx; -ry; -rz];

plot3(c1(1),c1(2),c1(3),'*r');
plot3(c2(1),c2(2),c2(3),'*r');
plot3(c3(1),c3(2),c3(3),'*r');
plot3(c4(1),c4(2),c4(3),'*r');
plot3(c5(1),c5(2),c5(3),'*r');
plot3(c6(1),c6(2),c6(3),'*r');
plot3(c7(1),c7(2),c7(3),'*r');
plot3(c8(1),c8(2),c8(3),'*r');