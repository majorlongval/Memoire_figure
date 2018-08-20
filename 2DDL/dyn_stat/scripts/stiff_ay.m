%% Script to determine best ay parameters for average cx and cy
clear all; close all; clc;
addpath('functions');

%% Arbitrary parameters giving the "size of the mechanism"
L = 5;
l = 1;


% The height at which the test is done :
x = L;

ay = linspace(-l,l,51);

cy = 0;
cx = 0;
for i = 1:length(ay)
            geo = pack_geo(ay(i),cy,cx,L,l);
            ymin = y_min_stat(geo);
            stat_WS0(i) = L-ymin;
            if ~(stat_WS0(i) == NaN);
                yv = linspace(ymin,L,50);
                for m =1:length(yv)
                    RLv(m) = rigid_lat([x,yv(m)],geo);
                    Rrv(m) = rigid_rot([x,yv(m)],geo);
                end
                RL(i) = mean(RLv);
                RR(i) = mean(Rrv);
            else
                RL(i) = NaN;
                RR(i) = NaN;
            end
    RL_ay_mean0(i) = RL(i);
    
    RR_ay_mean0(i) = RR(i);
end

 
cy = linspace(-l/4,l/4,51);
cx = cy;


RLv = zeros(1,51);
RL  = zeros(51,51,11);
RR  = zeros(51,51,11);
for i = 1:length(ay)
    for j = 1:length(cy)
        for k = 1:length(cx)
            geo = pack_geo(ay(i),cy(j),cx(k),L,l);
            ymin = y_min_stat(geo);
            stat_WS1(k,j,i) = L-ymin;
            if ~(stat_WS1(k,j,i) == NaN);
                yv = linspace(ymin,L,50);
                for m =1:length(yv)
                    RLv(m) = rigid_lat([x,yv(m)],geo);
                    RRv(m) = rigid_rot([x,yv(m)],geo);
                end
                RL(j,k,i) = mean(RLv);
                RR(j,k,i) = mean(RRv);
            else
                RL(j,k,i) = NaN;
                RR(j,k,i) = NaN;
            end
        end
    end
    RL_ay = RL(:,:,i);
    RL_ay_mean1(i) = mean(RL_ay(:));
    
    RR_ay = RR(:,:,i);
    RR_ay_mean1(i) = mean(RR_ay(:));
    
    stat_WS1_i = stat_WS1(:,:,i);
    stat_WS1_m(i) = mean(stat_WS1_i(:));
end

% cy = linspace(-1*l/2,l/2,51);
% cx = cy;
% 
% 
% RLv = zeros(1,51);
% RL  = zeros(51,51,11);
% RR  = zeros(51,51,11);
% for i = 1:length(ay)
%     for j = 1:length(cy)
%         for k = 1:length(cx)
%             geo = pack_geo(ay(i),cy(j),cx(k),L,l);
%             ymin = y_min_stat(geo);
%             stat_WS2(k,j,i) = L-ymin;
%             if ~(stat_WS2(k,j,i) == NaN);
%                 yv = linspace(ymin,L,50);
%                 for m =1:length(yv)
%                     RLv(m) = rigid_lat([x,yv(m)],geo);
%                     RRv(m) = rigid_rot([x,yv(m)],geo);
%                 end
%                 RL(j,k,i) = mean(RLv);
%                 RR(j,k,i) = mean(RRv);
%             else
%                 RL(j,k,i) = NaN;
%                 RR(j,k,i) = NaN;
%             end
%         end
%     end
%     RL_ay = RL(:,:,i);
%     RL_ay_mean2(i) = mean(RL_ay(:));
%     
%     RR_ay = RR(:,:,i);
%     RR_ay_mean2(i) = mean(RR_ay(:));
%     
%     stat_WS2_i = stat_WS2(:,:,i);
%     stat_WS2_m(i) = mean(stat_WS2_i(:));
% end    
%     
%     
% cy = linspace(-3*l/4,3*l/4,51);
% cx = cy;
% 
% 
% RLv = zeros(1,51);
% RL  = zeros(51,51,11);
% RR  = zeros(51,51,11);
% for i = 1:length(ay)
%     for j = 1:length(cy)
%         for k = 1:length(cx)
%             geo = pack_geo(ay(i),cy(j),cx(k),L,l);
%             ymin = y_min_stat(geo);
%             stat_WS3(k,j,i) = L-ymin;
%             if ~(stat_WS3(k,j,i) == NaN);
%                 yv = linspace(ymin,L,50);
%                 for m =1:length(yv)
%                     RLv(m) = rigid_lat([x,yv(m)],geo);
%                     RRv(m) = rigid_rot([x,yv(m)],geo);
%                 end
%                 RL(j,k,i) = mean(RLv);
%                 RR(j,k,i) = mean(RRv);
%             else
%                 RL(j,k,i) = NaN;
%                 RR(j,k,i) = NaN;
%             end
%         end
%     end
%     RL_ay = RL(:,:,i);
%     RL_ay_mean3(i) = mean(RL_ay(:));
%     
%     RR_ay = RR(:,:,i);
%     RR_ay_mean3(i) = mean(RR_ay(:));
%     
%         stat_WS3_i = stat_WS3(:,:,i);
%     stat_WS3_m(i) = mean(stat_WS3_i(:));
% end
% cy = linspace(-7*l/8,7*l/8,51);
% cx = cy;
% 
% 
% RLv = zeros(1,51);
% RL  = zeros(51,51,11);
% RR  = zeros(51,51,11);
% for i = 1:length(ay)
%     for j = 1:length(cy)
%         for k = 1:length(cx)
%             geo = pack_geo(ay(i),cy(j),cx(k),L,l);
%             ymin = y_min_stat(geo);
%             stat_WS4(k,j,i) = L-ymin;
%             if ~(stat_WS4(k,j,i) == NaN);
%                 yv = linspace(ymin,L,50);
%                 for m =1:length(yv)
%                     RLv(m) = rigid_lat([x,yv(m)],geo);
%                     RRv(m) = rigid_rot([x,yv(m)],geo);
%                 end
%                 RL(j,k,i) = mean(RLv);
%                 RR(j,k,i) = mean(RRv);
%             else
%                 RL(j,k,i) = NaN;
%                 RR(j,k,i) = NaN;
%             end
%         end
%     end
%     RL_ay = RL(:,:,i);
%     RL_ay_mean4(i) = mean(RL_ay(:));
%     
%     RR_ay = RR(:,:,i);
%     RR_ay_mean4(i) = mean(RR_ay(:));
%     
%     stat_WS4_i = stat_WS4(:,:,i);
%     stat_WS4_m(i) = mean(stat_WS4_i(:));
% end
% 
% % Finding the max
% [max_RL_ay_mean0,i_max_RL_ay_mean0] = max(RL_ay_mean0);
% [max_RL_ay_mean1,i_max_RL_ay_mean1] = max(RL_ay_mean1);
% [max_RL_ay_mean2,i_max_RL_ay_mean2] = max(RL_ay_mean2);
% [max_RL_ay_mean3,i_max_RL_ay_mean3] = max(RL_ay_mean3);
% [max_RL_ay_mean4,i_max_RL_ay_mean4] = max(RL_ay_mean4);
% 
% [max_RR_ay_mean0,i_max_RR_ay_mean0] = max(RR_ay_mean0);
% [max_RR_ay_mean1,i_max_RR_ay_mean1] = max(RR_ay_mean1);
% [max_RR_ay_mean2,i_max_RR_ay_mean2] = max(RR_ay_mean2);
% [max_RR_ay_mean3,i_max_RR_ay_mean3] = max(RR_ay_mean3);
% [max_RR_ay_mean4,i_max_RR_ay_mean4] = max(RR_ay_mean4);
% 
% [max_WS0,i_max_WS0] = max(stat_WS0);
% [max_WS1,i_max_WS1] = max(stat_WS1_m);
% [max_WS2,i_max_WS2] = max(stat_WS2_m);
% [max_WS3,i_max_WS3] = max(stat_WS3_m);
% [max_WS4,i_max_WS4] = max(stat_WS4_m);
% 
% [max_prod0,i_max_prod0] = max(RR_ay_mean0.*RL_ay_mean0.*stat_WS0);
% [max_prod1,i_max_prod1] = max(RR_ay_mean1.*RL_ay_mean1.*stat_WS1_m);
% [max_prod2,i_max_prod2] = max(RR_ay_mean2.*RL_ay_mean2.*stat_WS2_m);
% [max_prod3,i_max_prod3] = max(RR_ay_mean3.*RL_ay_mean3.*stat_WS3_m);
% [max_prod4,i_max_prod4] = max(RR_ay_mean4.*RL_ay_mean4.*stat_WS4_m);
% 
% 
% figure
% plot(ay,RL_ay_mean0/9.81,'-m');
% hold on;
% plot(ay,RL_ay_mean1/9.81,'-b');
% plot(ay,RL_ay_mean2/9.81,'-r');
% plot(ay,RL_ay_mean3/9.81,'-k');
% plot(ay,RL_ay_mean4/9.81,'-g');
% plot(ay(i_max_RL_ay_mean0),max_RL_ay_mean0/9.81,'*k');
% plot(ay(i_max_RL_ay_mean1),max_RL_ay_mean1/9.81,'*k');
% plot(ay(i_max_RL_ay_mean2),max_RL_ay_mean2/9.81,'*k');
% plot(ay(i_max_RL_ay_mean3),max_RL_ay_mean3/9.81,'*k');
% plot(ay(i_max_RL_ay_mean4),max_RL_ay_mean4/9.81,'*k');
% grid minor;
% xlabel('ay');
% ylabel('tauyg');
% title('title');
% legend('l0','l4','l2','3l4','7l8');
% 
% 
% figure
% plot(ay,RR_ay_mean0/9.81,'-m');
% hold on;
% plot(ay,RR_ay_mean1/9.81,'-b');
% plot(ay,RR_ay_mean2/9.81,'-r');
% plot(ay,RR_ay_mean3/9.81,'-k');
% plot(ay,RR_ay_mean4/9.81,'-g');
% plot(ay(i_max_RR_ay_mean0),max_RR_ay_mean0/9.81,'*k');
% plot(ay(i_max_RR_ay_mean1),max_RR_ay_mean1/9.81,'*k');
% plot(ay(i_max_RR_ay_mean2),max_RR_ay_mean2/9.81,'*k');
% plot(ay(i_max_RR_ay_mean3),max_RR_ay_mean3/9.81,'*k');
% plot(ay(i_max_RR_ay_mean4),max_RR_ay_mean4/9.81,'*k');
% grid minor;
% xlabel('ay');
% ylabel('taurg');
% title('title');
% legend('l0','l4','l2','3l4','7l8');
% 
% figure
% plot(ay,(RR_ay_mean0.*stat_WS0.*RL_ay_mean0),'-m');
% hold on;
% plot(ay,(RR_ay_mean1.*stat_WS1_m.*RL_ay_mean1),'-b');
% plot(ay,(RR_ay_mean2.*stat_WS2_m.*RL_ay_mean2),'-r');
% plot(ay,(RR_ay_mean3.*stat_WS3_m.*RL_ay_mean3),'-k');
% plot(ay,(RR_ay_mean4.*stat_WS4_m.*RL_ay_mean4),'-g');
% plot(ay(i_max_prod0),max_prod0,'*k');
% plot(ay(i_max_prod1),max_prod1,'*k');
% plot(ay(i_max_prod2),max_prod2,'*k');
% plot(ay(i_max_prod3),max_prod3,'*k');
% plot(ay(i_max_prod4),max_prod4,'*k');
% grid minor;
% xlabel('ay');
% ylabel('taurg');
% title('title');
% legend('l0','l4','l2','3l4','7l8');
% 
% figure
% plot(ay,stat_WS0,'-m');
% hold on;
% plot(ay,stat_WS1_m,'-b');
% plot(ay,stat_WS2_m,'-r');
% plot(ay,stat_WS3_m,'-k');
% plot(ay,stat_WS4_m,'-g');
% grid minor;
% xlabel('ay');
% ylabel('taurg');
% title('title');
% legend('l0','l4','l2','3l4','7l8');
% 
