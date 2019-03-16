addpath('C:\Users\jorda\OneDrive\Desktop\Memoire\Memoire_figure\3DDL\functions');

alpha  = pi/6;     % rad
R      = 0.7033;   % m 
r      = 0.2;      % m
g      = 9.81;     % m/s^2

phic = 0;
m_cell = 0.174;    % kg mass of cellphone
m_eff  = 0.523;    % kg mass of end effector
m_w    = 0.376;    % mass of extra weight
m  = m_cell + m_eff + m_w;

h_cell = 0.035;    % m
h_w    = 0.015;    % m
h_eff  = 0.04;     % m 

%% Calculating rc 
rc = m_w*r/m;

%% Calculating hc
hc = (m_cell*h_cell + h_w*m_w + h_eff*m_eff)/m;


% Setting the heights
z1 = 2; 
z2 = 2.5;

% Setting the wc
% Range of values for the external wrench and the mass
Mx = 0;     % Nm
My = 0;      % Nm
Mz = 0;    % Nm
fx = 0;     % N
fy = 9.81*0.138;      % N
fz = 0;          % N


% Defining all combinations of wrenches ( can define multiple fo mult.
% plot)
wc = wrench_comb(fx,fy,fz,Mz,My,Mz);

pii_0 = WFW_at_h(alpha, R, r, rc, phic, hc,m,zeros(1,6),z1);
pii_1 = WFW_at_h(alpha, R, r, rc, phic, hc,m,wc,z1);
pii_2 = WFW_at_h(alpha, R, r, rc, phic, hc,m,wc,z2);

[K0,vol0] = convhull(pii_0(:,1),pii_0(:,2));
[K1,vol1] = convhull(pii_1(:,1),pii_1(:,2));
[K2,vol2] = convhull(pii_2(:,1),pii_2(:,2));


% Plotting the points 
% figure;
% grid on;
% hold;

% plot3(pii_0(K0,1),pii_0(K0,2),pii_0(K0,3),'-g');
% plot3(pii_1(K1,1),pii_1(K1,2),pii_1(K1,3),'-b');
% plot3(pii_2(K2,1),pii_2(K2,2),pii_2(K2,3),'-r');
% plot3(bc1(1),bc1(2),bc1(3),'*b');
% plot3(bc2(1),bc2(2),bc2(3),'*r');
% plot3(bc0(1),bc0(2),bc0(3),'*g');


Plotting the other points of the trajectory
n1 = length(pii_1(:,1));
bc1 = mean(pii_1);
bc2 = mean(pii_2);
bc0 = mean(pii_0);
eta = 0.6:0.1:0.9;
points_0 = [];
points_1 = [];
points_2 = [];
for i = 1:length(eta)
    points_0 = [points_0;bc0+eta(i)*(pii_0(K0,:)-bc0)];
    points_1 = [points_1;bc1+eta(i)*(pii_1(K1,:)-bc1)];
    points_2 = [points_2;bc2+eta(i)*(pii_2(K2,:)-bc2)];
end
points_0 = [bc0;points_0];
points_1 = [bc1;points_1];
points_2 = [bc2;points_2];
points  = [points_0;points_1;points_2];
setup_points = [points(:,1),-1*points(:,2),points(:,3)-0.3170];

plot3(points_0(:,1),points_0(:,2),points_0(:,3),'--g');
plot3(points_1(:,1),points_1(:,2),points_1(:,3),'--b');
plot3(points_2(:,1),points_2(:,2),points_2(:,3),'--r');

figure;
view(3);
grid on;
hold on;
plot3(points(:,1),points(:,2),points(:,3),'--k');

% save points_funkytown.mat setup_points -v4
function wc = wrench_comb(fx,fy,fz,Mx,My,Mz)
wc = [];
for a = 1:length(fx)
    for b =1:length(fy)
        for c = 1:length(fz)
            for d = 1:length(Mx)
                for e = 1:length(My)
                    for f = 1:length(Mz)
                        wc = [wc;...
                            fx(a),fy(b),fz(c),Mx(d),My(e),Mz(f)];
                    end
                end
            end
        end
    end
end
end