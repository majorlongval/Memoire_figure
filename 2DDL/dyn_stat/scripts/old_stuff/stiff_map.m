%% This script maps the latteral and rotational stifness of the mechanism
clear all; close all; clc;
addpath('functions');

% Setting the geometric parameters
L = 5;
l = 1;
cx = 0;

xpos = 2;  %linspace(1,5,100);
%ypos = linspace(-L, L, 100);

aya = linspace(-l,l,50);
cya = linspace(-l,l,50);


for i =1:length(aya)
    B2   = aya(i)-l;
    B3   = -(aya(i)+l);
    for j =1:length(cya)
        A3 = aya(i)*(L-cya(j))-2*L*cya(j)-l*(aya(i)+L);
        A2 = -aya(i)*(L-cya(j))+2*L*cya(j)-l*(aya(i)+L);
        gp = struct('L',L,'l',l,'ay',aya(i),'cy',cya(j),'cx',cx);
        ymin  = max(-A3/B3,-A2/B2);
        S(i,j) = L - max(-A3/B3,-A2/B2);
        ypos = linspace(ymin,L,100);
        for k = 1:length(ypos)
            [stiff_lat(k),stiff_rot(k)] = Rigid([xpos,ypos(k)],gp);
        end
        stiff_lat_m(i,j) = mean(stiff_lat);
        stiff_rot_m(i,j) = mean(stiff_rot);
    end
end

%  Trouver les max de chaque stiff et de S

[max_stiff_lat,i_max_stiff_lat]=max(stiff_lat_m(:));
[I_row_max_stiff_lat, I_column_max_stiff_lat] = ...
    ind2sub(size(stiff_lat_m),i_max_stiff_lat);

[max_stiff_rot,i_max_stiff_rot]=max(stiff_rot_m(:));
[I_row_max_stiff_rot, I_column_max_stiff_rot] = ...
    ind2sub(size(stiff_rot_m),i_max_stiff_rot);

[max_S,i_max_S]=max(S(:));
[I_row_max_S, I_column_max_S] = ...
    ind2sub(size(S),i_max_S);

max_vect_S = [];
for i =1:length(S(:,1))
    for j =1:length(S(1,:))
        if S(i,j) >= 2*L
           max_vect_S=[max_vect_S;i,j,S(i,j)];
        end
    end
end

norm_stiff_lat_m = stiff_lat_m/max_stiff_lat;
norm_stiff_rot_m = stiff_rot_m/max_stiff_rot;
norm_S = S/max_S;

norm_prod = norm_stiff_lat_m.*norm_stiff_rot_m.*norm_S;

prod = stiff_lat_m.*stiff_rot_m.*S;

[max_norm_prod,i_max_norm_prod]=max(norm_prod(:));
[I_row_norm_prod, I_column_norm_prod] = ...
    ind2sub(size(norm_prod),i_max_norm_prod);

[max_prod,i_max_prod]=max(prod(:));
[I_row_prod, I_column_prod] = ...
    ind2sub(size(prod),i_max_prod);

figure
C = contourf(aya,cya, norm_prod/max_norm_prod,20);
h = colorbar;
ylabel(h, 'prod')
xlabel('ay'); ylabel('cy');
grid minor;
hold on;
plot(aya(I_column_norm_prod),cya(I_row_norm_prod),'*r');


figure
contourf(aya,cya, norm_stiff_lat_m,0:0.1:1);
h = colorbar;
ylabel(h, 'stiff_lat')
xlabel('ay'); ylabel('cy');
grid minor;

figure
contourf(aya,cya, norm_stiff_rot_m,0:0.1:1);
h = colorbar;
ylabel(h, 'stiff_rot')
xlabel('ay'); ylabel('cy');
grid minor;

f4 = figure;
ax = gca;
ax.XTickMode = 'manual';
v = [0:0.1:1,1,1.05];
[C,h1] = contourf(aya,cya, S/(2*L),v);
% h = colorbar('Ticks',[min(S(:))/(2*L),0.1:0.1:1,max(S(:))/(2*L)]);
% ylabel(h, 'S')
xlabel('ay'); ylabel('cy');
xticks([-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1]);
grid on;
hold on;
p4 = plot([-1 1],[-1 1],'-r','LineWidth',2);
% clabel(C,h1,'manual');
title('Ltitle');
annotation('textarrow',[(-0.35+1)/2 (-0.55+1)/2],[(-0.7+1)/2 (-0.6+1)/2],'String','annot');


% print(f4,'/home/jordan/Documents/Maitrise/recherche/memoire_figures/2DDL/dyn_stat/figs/espace_stat_max','-depsc');

print(f4,'-depsc','-painters','/home/jordan/Documents/Maitrise/recherche/memoire_figures/2DDL/dyn_stat/figs/espace_stat_max.eps')
