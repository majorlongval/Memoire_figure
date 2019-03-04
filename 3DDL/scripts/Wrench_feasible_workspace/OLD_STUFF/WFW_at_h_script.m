clear all; close all;clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');


% Setting the arbitrary values
R     = 0.75;                 % m
r     = 0.2;                  % m
alpha = pi/6;                 % rad
rc    = 0.01;              % m
phic  = 0;                    % rad
hc    = 0;                    % m


% Setting the minimum and maximum heigth
zmin = 2;                     % m
zmax = 10;                    % m

% Range of values for the external wrench and the mass
Mx = [-1 1];                       % Nm
My = [-1 1];                       % Nm
Mz = [-1 1];                       % Nm
fx = [-1 1];                  % N
fy = [-1 1];                       % N
fz = 0;                       % N
m  = 5;                     % kg

% Defining all combinations of wrenches ( can define multiple fo mult.
% plot)
wc = wrench_comb(fx,fy,fz,Mz,My,Mz);
% R     = 0.75;                 % m
% r     = 0.2;                  % m
% alpha = pi/6;                 % rad
% rc    = 0.01;              % m
% phic  = 0;                    % rad
% hc    = 0;  % R     = 0.75;              % m
% % r     = 0.2;               % m
% % alpha = pi/6;              % rad
% % rc    = 0.00001;           % m
% % phic  = 0;                 % rad
% % hc    = 0.00001;            % m
% 
% 
% z = 2;                     % m
% % Range of values for the external wrench and the mass
% Mx = [-1 1];                       % Nm
% My = [-1 1];                       % Nm
% Mz = [-1 1];                       % Nm
% fx = [-1 1];                  % N
% fy = [-1 1];% Mx = [-0.2*.376*9.81*cos(pi/4),0.2*.376*9.81*cos(pi/4)];                    % Nm
% % My = [-0.2*.376*9.81*cos(pi/4),0.2*.376*9.81*cos(pi/4)];          % Nm
% % Mz = 0;                    % Nm
% % fx = 0;                    % N
% % fy = 0;                    % N
% fz =0;% fz = .376*9.81;             % N
% m = 10;% m  = 0.532+0.1741;                    % kg
% 
% wc = wrench_comb(fx,fy,fz,Mx,My,Mz);

corners = WFW_at_h(alpha,R,r,rc,phic,hc,m,wc,zmin);
phys_corners = [corners(:,1),-1*corners(:,2),ones(length(corners(:,1)),1)*1.683]

bary_phys_corners = mean(phys_corners);


plot(phys_corners(:,1),phys_corners(:,2),'*r');
phys_corners_05 = [0.5*phys_corners(:,1:2),phys_corners(:,3)];



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