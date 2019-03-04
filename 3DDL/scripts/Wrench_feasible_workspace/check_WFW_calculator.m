% Script to check if the WFW is good.
clear all;close all; clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

% Parameters
R     = 1;         % m
r     = 0.2;       % m
alpha = pi/6;         % rad
rc    = 0;         % m
pc  = 0;      % rad
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





[corners,volume,K] = WFW_Calculator(R,r,alpha,rc,pc,hc,m,wrench_comb,zr);
for i = 1:length(corners(:,1))
    planes_of_corners(i,:) = sort(corners(i,1:3));
end
for i = 1:length(planes_of_corners(:,1))
    test = planes_of_corners;
    test(i,:) = [];
    if ismember(planes_of_corners(i,:),test,'rows')
        corners(i,:) = NaN*ones(1,length(corners(1,:)));
        planes_of_corners(i,:) = NaN*ones(1,length(planes_of_corners(1,:)));
    end
end
corners(any(isnan(corners),2),:) = [];
corners = sortrows(corners,6);


T = 200;
t = linspace(0,T,100);
tau = t/T;
pos = [];
acc = [];
tt = [];
q1 = 6*tau.^5-15*tau.^4+10*tau.^3;
q3 = 120*tau.^3-180*tau.^2+60*tau;

for i = 1:length(corners(:,1))-1
    tt = [tt,i*t(end)+t];
    pos = [pos,corners(i,4:6)'+(corners(i+1,4:6)-corners(i,4:6))'*q1];
    acc = [acc,(corners(i+1,4:6)-corners(i,4:6))'*q3/(T^2)];
end


plot3(corners(:,4),corners(:,5),corners(:,6),'.r');
hold on; grid on;
plot3(pos(1,:),pos(2,:),pos(3,:),'-b');



% figure;
% plot(tt,pos(1,:));
% hold on;
% plot(tt,pos(2,:));
% plot(tt,pos(3,:));
% 
% figure;
% plot(tt,acc(1,:));
% hold on;
% plot(tt,acc(2,:));
% plot(tt,acc(3,:));
% 
% 
% 
%%%%%%%%%%%%% CHECHING IF THE TENSION IS POSITIVE %%%%%%%%%%%%%%%%%%%%
g = 9.81;
a11 = r*[0; 1; 0];
a12 = -a11;
ct = cos(2*pi/3); st = sin(2*pi/3);
Qs = [ct -st 0;st ct 0;0 0 1];
a21 = Qs*a11;
a22 = -a21;
a31 = Qs*a21;
a32 = -a31;
alpha = pi/6;
R1 = R*[cos(alpha);sin(alpha);0];
R2 = Qs*R1;
R3 = Qs*R2;






b11 = R1+a11;
b12 = R1+a12;
b21 = R2+a21;
b22 = R2+a22;
b31 = R3+a31;
b32 = R3+a32;

for j =1:length(wrench_comb(:,1))
    wc = wrench_comb(j,:);
    for i =1:length(pos)
        rho1(:,i) = pos(:,i)-R1;
        rho2(:,i) = rho1(:,i);
        rho3(:,i) = pos(:,i)-R2;
        rho4(:,i) = rho3(:,i);
        rho5(:,i) = pos(:,i)-R3;
        rho6(:,i) = rho5(:,i);
        e1(:,i) = rho1(:,i)/norm(rho1(:,i));
        e2(:,i) = rho3(:,i)/norm(rho3(:,i));
        e3(:,i) = rho5(:,i)/norm(rho5(:,i));
        delta11 = cross(a11-c,e1(:,i));
        delta12 = cross(a12-c,e1(:,i));
        delta21 = cross(a21-c,e2(:,i));
        delta22 = cross(a22-c,e2(:,i));
        delta31 = cross(a31-c,e3(:,i));
        delta32 = cross(a32-c,e3(:,i));
        
        M(:,:,i) = [rho1(:,i)/norm(rho1(:,i)) rho2(:,i)/norm(rho2(:,i)) rho3(:,i)/norm(rho3(:,i))...
            rho4(:,i)/norm(rho4(:,i)) rho5(:,i)/norm(rho5(:,i)) rho6(:,i)/norm(rho6(:,i));
            delta11 delta12 delta21 delta22 delta31 delta32];
        
        gamma(:,i) = [wc(1)-m*acc(1,i);wc(2)-m*acc(2,i);wc(3)-m*acc(3,i)+m*g;wc(4);wc(5);wc(6)];
        tens(:,i,j) = M(:,:,i)\gamma(:,i);
    end
end

for j =1:length(wrench_comb(:,1))
    figure
    for i =1:6
        plot(tt,tens(i,:,j));
        hold on;
    end
end