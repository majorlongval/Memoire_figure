% Script pour montrer si l'ensemble des torseurs à l'effecteur est à
% l'intérieur de l'enveloppe ouverte des torseurs disponibles.
clear all; close all; clc;

% Position to be tested
x = 0; y = 0; z = 5;

% Geometric disposition
R = 1;
r = 0.2;
alpha = pi/4;
rc = 0.04; phic = pi/8; hc = 0.1;

m = 1;
g = 9.81;
% Geometric terms
[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, rc, phic, hc);

% Defining the planes
px = ay*y+az*z+ak;
py = bx*x+bz*z+bk;
pz = cx*x+cy*y+ck;

% Generating the data and plotting the planes
tx = [1 -1 -1 1]; % Generate data for x vertices
ty = [1 1 -1 -1]; % Generate data for y vertices
for i =1:length(px)
    tz(:,:,i) = (-1/pz(i))*(px(i)*tx + py(i)*ty + m*g*pz(i));
    ph= patch(tx,ty,tz(:,:,i),[0.5,0.9,0.5]);
    hold on;
    ph.FaceAlpha = 0.4;
    ph.LineWidth = 1.5;
end
view(3);
grid on;
xlim([-1 1]);
ylim([-1 1]);


% Plotting the intersection lines
line_comb = combnk(1:6,2);
z_line = linspace(-40,40,100);
for i = 1:length(line_comb(:,1))
   for j =1:length(z_line)
      txy(:,j,i) =  [px(line_comb(i,1)),py(line_comb(i,1)),pz(line_comb(i,1));...
                     px(line_comb(i,2)),py(line_comb(i,2)),pz(line_comb(i,2));...
                     0,                 0,                 1]\...
                     [-m*g*pz(line_comb(i,1));-m*g*pz(line_comb(i,2));-z_line(j)-m*g];
   end
   plot3(txy(1,:,i),txy(2,:,i),txy(3,:,i),'-k','LineWidth',1.5);
   hold on;
end

% Plotting the required wrench set
txmin = -0.5; txmax = 0.5;
tymin = -0.5; tymax = 0.5;
tzmin = -8; tzmax = 12;

ph1 = patch([txmax txmax txmax txmax],...
      [tymin tymax tymax tymin],...
      [tzmin tzmin tzmax tzmax],[0.2,0.2,0.9]);
ph1.FaceAlpha = 0.5;
ph1.LineWidth = 1.5;

ph2 = patch([txmin txmax txmax txmin],...
      [tymin tymin tymin tymin],...
      [tzmin tzmin tzmax tzmax],[0.2,0.2,0.9]);
ph2.FaceAlpha = 0.5;
ph2.LineWidth = 1.5;

ph3 = patch([txmin txmin txmin txmin],...
      [tymin tymax tymax tymin],...
      [tzmax tzmax tzmin tzmin],[0.2,0.2,0.9]);
ph3.FaceAlpha = 0.5;
ph3.LineWidth = 1.5;

ph4 = patch([txmin txmax txmax txmin],...
      [tymax tymax tymax tymax],...
      [tzmax tzmax tzmin tzmin],[0.2,0.2,0.9]);
ph4.FaceAlpha = 0.5;
ph4.LineWidth = 1.5;

ph5 = patch([txmax txmin txmin txmax],...
      [tymin tymin tymax tymax],...
      [tzmax tzmax tzmax tzmax],[0.2,0.2,0.9]);
ph5.FaceAlpha = 0.5;
ph5.LineWidth = 1.5;

ph6 = patch([txmax txmin txmin txmax],...
      [tymin tymin tymax tymax],...
      [tzmin tzmin tzmin tzmin],[0.2,0.2,0.9]);
ph6.FaceAlpha = 0.5;
ph6.LineWidth = 1.5;

xlabel('XLABEL');
ylabel('YLABEL');
zlabel('ZLABEL');

fig_h = gca;


print('bad_tors','-dsvg','-painters');