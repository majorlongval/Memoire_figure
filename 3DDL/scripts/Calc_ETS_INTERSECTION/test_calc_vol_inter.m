%%% TESTER la function Calc_vol_inter
clear all;clc;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
alpha = 0;
R = 1;
r = 0.2;
r_c = 0.010;
phi_c = 0;
h_c = 0;
tx = 0; 
ty = -2;
tz = 0;
m  = 3;
Mx = 2;
My = 0;
Mz = 0;
max_z = 10;

figure
discret_z = 20;
% Calcul de l'ETS INTER
[K_inter,vol_inter,seg_inter] = ...
calc_vol_INTER(alpha, R, r, r_c, phi_c,h_c, tx, ty, tz, m, Mx,My, Mz,max_z,discret_z);
plot3(seg_inter(:,1),seg_inter(:,2),seg_inter(:,3));
tri_handle_INTER = ...
    trisurf(K_inter,seg_inter(:,1),seg_inter(:,2),seg_inter(:,3),...
                             ones(length(seg_inter(:,1)),3)*[1;0;0]);
tri_handle_INTER.FaceAlpha = 0.8;
% 
% Calcul de l'ETST
% [K_ETST,vol_ETST,seg_ETST] = ...
% calc_vol_ETST(alpha, R, r, r_c, phi_c,h_c, tx, ty, tz, m, Mx,My, Mz,max_z,discret_z);
% % 
% plot3(seg_ETST(:,1),seg_ETST(:,2),seg_ETST(:,3),'-r');
% hold on;
% handle_red = plot3(seg_ETST(:,1),seg_ETST(:,2),seg_ETST(:,3),'*r');
% % hold on;
% % tri_handle_ETST = ...
% %     trisurf(K_ETST,seg_ETST(:,1),seg_ETST(:,2),seg_ETST(:,3),...
% %                              ones(length(seg_ETST(:,1)),3)*[1;0;0]);
% % tri_handle_ETST.FaceAlpha = 0.0;
% % tri_handle_ETST.EdgeColor = [1,0,0];
% % tri_handle_ETST.LineWidth = 2;
% % % Calcul de l'ETS
% [K_ETS,vol_ETS,seg_ETS] = ...
% calc_vol_ETST(alpha, R, r, r_c, phi_c,h_c, 0, 0, 0, m, 0,0, 0,max_z,discret_z);
% plot3(seg_ETS(:,1),seg_ETS(:,2),seg_ETS(:,3),'-b');
% hold on;
% handle_blue  = plot3(seg_ETS(:,1),seg_ETS(:,2),seg_ETS(:,3),'*b');
% tri_handle_ETS = ...
%     trisurf(K_ETS,seg_ETS(:,1),seg_ETS(:,2),seg_ETS(:,3),...
%                              ones(length(seg_ETS(:,1)),3)*[0.2;0.2;0.2]);
% tri_handle_ETS.FaceAlpha = 0.0;
% tri_handle_ETS.EdgeColor = [0,1,0];
% tri_handle_ETS.LineWidth = 2;
% % 
% % 
% % Calcul de l'ETS  max 
% [K_ETS_max,vol_ETS_max,seg_ETS_max] = ...
% calc_vol_ETST(alpha, R, r, 0, 0,0, 0, 0, 0, m, 0,0, 0,max_z,discret_z);
% tri_handle_ETS_max = ...
%     trisurf(K_ETS_max,seg_ETS_max(:,1),seg_ETS_max(:,2),seg_ETS_max(:,3),...
%                              ones(length(seg_ETS_max(:,1)),3)*[0.9;0.9;0.9]);
% tri_handle_ETS_max.FaceAlpha = 0.0;
% tri_handle_ETS_max.EdgeColor = [0,0,1];
% tri_handle_ETS_max.LineWidth = 2;

%  % Tracer les espaces
%  figure;
%  zr = linspace(0,max_z,discret_z);
% % for i =1:length(seg_inter)
% %     seg_int = cell2mat(seg_inter(i));
% %     plot3([seg_int(:,1);seg_int(1,1)],...
% %           [seg_int(:,2);seg_int(1,2)],...
% %           ones(length(seg_int(:,1))+1,1)*zr(i),'-r');
% %       hold on;
% % end
% for i =1:length(seg_ETST)
%     seg_ETST_T = cell2mat(seg_ETST(i));
%     plot3([seg_ETST_T(:,1);seg_ETST_T(1,1)],...
%           [seg_ETST_T(:,2);seg_ETST_T(1,2)],...
%           ones(length(seg_ETST_T(:,1))+1,1)*zr(i),'-b');
%       hold on;
% end
% for i =1:length(seg_ETS)
%     seg_ETS_T = cell2mat(seg_ETS(i));
%     plot3([seg_ETS_T(:,1);seg_ETS_T(1,1)],...
%           [seg_ETS_T(:,2);seg_ETS_T(1,2)],...
%           ones(length(seg_ETS_T(:,1))+1,1)*zr(i),'-g');
%       hold on;
% end
% for i =1:length(seg_ETS_max)
%     seg_ETS_max_T = cell2mat(seg_ETS_max(i));
%     plot3([seg_ETS_max_T(:,1);seg_ETS_max_T(1,1)],...
%           [seg_ETS_max_T(:,2);seg_ETS_max_T(1,2)],...
%           ones(length(seg_ETS_max_T(:,1))+1,1)*zr(i),'-y');
%       hold on;
% end

set(gca,'Zdir','reverse');
set(gca,'Xdir','reverse');
grid on;
xlabel('XXX');
ylabel('YYY');
zlabel('ZZZ');
title({'TITLE1','TITLE2'});
view([-79,34]);
%legend([handle_red,handle_blue],{'ETST','ETS'});
%saveas(gca,'Exemple_ETST.eps','epsc');
%[s_inter1,s_inter_2,s_inter3] = size(seg_inter);
% zr = linspace(0,max_z,discret_z);
% for i = 1:s_inter3
%     plot3([seg_inter(:,1,i);seg_inter(1,1,i)],...
%          [seg_inter(:,2,i);seg_inter(1,2,i)],...
%           ones(length(seg_inter(:,1,i))+1,1)*zr(i),'-r');
%      hold on;
% end