%% Calculating the volume and vertices of the wrench feasible workspace
% Author : Jordan M. Longval
% This is the third version of this script. I think it's the best so far.
clear all; close all; clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

% Setting the arbitrary values
R     = 1;         % m
r     = 0.2;       % m
alpha = 0;         % rad
rc    = 0.04;         % m
phic  = pi/6;      % rad
hc    = 0.5;         % m


% Setting the minimum and maximum heigth
zmin = 1;          % m
zmax = 10;         % m

% Range of values for the external wrench and the mass
Mx = [-2 2];     % Nm
My = [2 6];      % Nm
Mz = [1 1.5];    % Nm
fx = [-6 2];     % N
fy = [4 8];      % N
fz = 0;          % N
m  = 30;         % kg


% Defining all combinations of wrenches ( can define multiple fo mult.
% plot)
wc = wrench_comb(fx,fy,fz,Mz,My,Mz);




% Calculating all the possible planes
planes = [];
for k =1:length(wc(:,1))
    planes = [planes;Calc_planes(alpha,R,r,rc,phic,hc,m,...
        wc(k,:))];
end
planes = [planes;0 0 -1,zmax;0 0 1 -zmin];

% Defining the list of all projection lines at the plane z =zmin
L_l_i = [];
for i =1:length(wc(:,1))
    L_l_i = [L_l_i;Calc_lines(alpha, R, r, rc, phic, ...
        hc,m,zmin,wc(i,:))];
end

% Determining the points on the z=zmin plane which are good.
tic
pii = [];
for i =1:length(L_l_i(:,1))-1
    for j=i+1:length(L_l_i(:,1))
        temp = Calc_intersection(L_l_i(i,:), L_l_i(j,:));
        bool = 1;
        for k =1:length(L_l_i(:,1))
            if k ~= i && k~= j
                verif = L_l_i(k,1)*temp(1)+L_l_i(k,2)*temp(2)+L_l_i(k,3) >= 0;
                bool = bool*verif;
            end
        end
        if bool ==1
            pii = [pii;i,j,194,temp,zmin];
            plot3(temp(1),temp(2),temp(3));
            hold on;
            pause;
            
        end
    end
end
toc



% Determining all_the points on the line and determining the shortest
% segment
ishort = 1;
points = cell(length(pii(:,1)),1);
new_p = [];
for i =1:length(pii(:,1))
    dist = [];
    for j =1:length(planes(:,1))
        if ~ismember(j,pii(i,1:3))
            x = calc_plane_intersection(planes(pii(i,1),:),...
                planes(pii(i,2),:),...
                planes(j,:));
            if ~isempty(x) && x(3)> pii(i,6)
                temp = points{i};
                temp = [temp,x];
                points{i} = temp;
                dist = [dist;j,dist_2_points(pii(i,4:6),x)];
            end
        end
    end
    dist = sortrows(dist,2,'ascend');
    short = dist(1,:);
    pii = [pii; pii(i,1:2),short(1),calc_plane_intersection(planes(pii(i,1),:),...
                                           planes(pii(i,2),:),...
                                           planes(short(1),:))'];
end
pii = sortrows(pii,6);
m_planes = uint8(pii(:,1:3));
m_points = pii(:,4:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    FUNCTIONS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i =1:length(pii(:,1))
   plot3(pii(i,4),pii(i,5),pii(i,6),'.r','markerSize',20);% 
   hold on;
end
grid on;
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
function dist = dist_2_points(p1,p2)
    dist = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2+(p1(3)-p2(3))^2);
end