%% Calculating the boundaries of the wrench feasible workspace
% Author : Jordan M. Longval

% Initialising 
clear; clc; close all;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
% Setting the arbitrary values 
R     = 1;      % m
r     = 0.2;    % m
alpha = 0;   % rad
rc    = 0;   % m
phic  = pi/6;   % rad
hc    = 0;    % m 


% Setting the maximum z 
zmax = 10;


% Maximal and minimal cable tension  
fmax = 50; % N
fmin = 0.2;   % N 

% Range of values for the external wrench and the mass
Mx = [-0.4, -0.2];%[-0.02 0.02]; % Nm
My = [-0.01 0.05]; % Nm
Mz = [-0.1 0.1]; % Nm
fx = [-1  1]; % N
fy = [-1 1]; % N
fz = 0;%[-10  10]; % N
m  = 3; %kg

% Span of the z axis incrementation
z_span = linspace(0,zmax,17);

% Limits of the task workspace
center = [ 0  0  4];
amp = 0.2;
corners = center + amp*[1 -1 -1;
                       -1 -1 -1;
                        -1 -1 1;
                        1 -1 1;
                         1 1 1;
                        -1 1 1;
                       -1 1 -1;
                        1 1 -1];
% Max  and min cable lengths
[rhomax, rhomin] = calc_rhomax_rhomin(R, alpha, corners);                 

[K,vertices,volume] = vertices_WFW(R,r,alpha,rc,phic,hc,fx,fy,fz,Mx,My,Mz,m,z_span,rhomax,rhomin,fmin,fmax);

figure;

plot3(vertices(:,1),vertices(:,2),vertices(:,3),'*r');
hold on;
plot3(corners(:,1),corners(:,2),corners(:,3),'*b');

set(gca,'Zdir','reverse');
set(gca,'Xdir','reverse');
grid on;
xlabel('XXX');
ylabel('YYY');
zlabel('ZZZ');
title({'TITLE1','TITLE2'});
view([-79,34]);

% ub = [(pi/2)-0.001,r/2,pi/3,r];
% lb = [0,-r/2,0,-r];
% x0 = (ub+lb)/2;
% 
% x_sol = opti_WFW(R, r,  ub, lb, x0,fx,fy,fz,Mx,My,Mz,m,z_span)

function x_sol = opti_WFW(R, r,  ub, lb, x0,fx,fy,fz,Mx,My,Mz,m,z_span)
        fun = @(x)-vertices_WFW(R,r,x(1),x(2),x(3),x(4),fx,fy,fz,Mx,My,Mz,m,z_span);
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        nonlcon = [];
        %options = optimoptions('fmincon','Algorithm','sqp');
        options = optimoptions('fmincon','Display','iter');
        x_sol = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end


function [K,vertices,volume] = vertices_WFW(R,r,alpha,rc,phic,hc,fx,fy,fz,Mx,My,Mz,m,z_span,rhomax,rhomin,fmin,fmax)
% Making all the possible combinaitions of extreme wrench values

vertices = [];
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

for z =1:length(z_span)
% Calculating all the lines 
List_lines = [];
for i =1:length(wrench_comb(:,1))
    List_lines = [List_lines;Calc_lines(alpha, R, r, rc, phic, ...
                                         hc,m,z_span(z),wrench_comb(i,:),rhomax,rhomin,fmin,fmax)];
end

% Calculating all the intersection points 
intersection = cell(length(List_lines(:,1)),1);
counter = 0;
for i =1:length(List_lines(:,1))-1
    for j =i+1:length(List_lines(:,1))
        p_inter = Calc_intersection(List_lines(i,:), List_lines(j,:));
        temp_i = intersection{i};
        temp_i = [temp_i;p_inter];
        intersection{i} = temp_i;
        temp_j = intersection{j};
        temp_j = [temp_j;p_inter];
        intersection{j} = temp_j;
        counter = counter+1;
    end
end

% Ordering the points according to a lexicographic order and checking which
% segments of each line respect all the conditions
roundn = @(x,n) round(x.*10.^n)./10.^n;
for i =1:length(intersection)
    temp = intersection{i};
    [Xp, Yp] = vec_acw_order_lin(temp(:,1),temp(:,2));
    for j =1:length(Xp)
        for k =1:length(List_lines(:,1))
            bool(k) = roundn(List_lines(k,1)*Xp(j)+...
                             List_lines(k,2)*Yp(j)+List_lines(k,3),5)>=0;
        end
        if all(bool== 1)
            vertices = [vertices;Xp(j) Yp(j) z_span(z)];
        end
    end
end
end
vertices = roundn(vertices,4);
vertices = unique(vertices,'rows');
vertices = sortrows(vertices,3);

[K,volume] = convhull(vertices);

end
        

function [rhomax, rhomin] = calc_rhomax_rhomin(R, alpha, corners)
Qs = [cos(2*pi/3) -sin(2*pi/3) 0;...
    sin(2*pi/3) cos(2*pi/3) 0;
    0 0 1];

RR(:,1) = R*[cos(alpha);sin(alpha);0];
RR(:,2) = RR(:,1);
RR(:,3) = Qs*RR(:,2);
RR(:,4) = Qs*RR(:,2);
RR(:,5) = Qs*RR(:,4);
RR(:,6) = Qs*RR(:,4);

rhomin = [1000;1000;1000;1000;1000;1000];
rhomax = [0; 0; 0; 0; 0; 0];
for i =1:length(rhomin)
    for j =1:length(corners(:,1))
        if norm(corners(j,:)'-RR(:,i))<=rhomin(i)
            rhomin(i) = norm(corners(j,:)'-RR(:,i));
        elseif norm(corners(j,:)'-RR(:,i))>=rhomax(i)
            rhomax(i) = norm(corners(j,:)'-RR(:,i));
        end
    end
end
end