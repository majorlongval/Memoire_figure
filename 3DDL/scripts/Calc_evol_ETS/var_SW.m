%% Script to analyse the influence of the parameters on the area of SW
clear all;close all;clc;
alpha =linspace(0,(pi/2)-0.001,50);
R = 1;
r = 0.5;
rc = linspace(0.001,r,10);
pc = linspace(0,pi/3,50);
alpha_m_mean = 0;
pc_m_mean = 0;
bool = 0;
while bool == 0
for k =1:length(rc)
    for i =1:length(alpha)
        for j =1:length(pc)
            area(j,i) = calc_area_SW_new(alpha(i),R,r,rc(k),pc(j));
            %             display('j= ');
            %             display(j);
            %             display('\n');
        end
        %         display('i= ');
        %             display(i);
        %             display('\n');
    end
    display(k);
    [M,I] = max(area(:))
    [row,col] = ind2sub(size(area),I);
    alpha_m(k) = alpha(col);
    pc_m(k) = pc(row);
%     figure
%     contour(alpha,pc,area);
%     hold on;
%     plot(alpha_m(k),pc_m(k),'*r');
    
end
display(alpha_m_mean - mean(alpha_m));
display(pc_m_mean - mean(pc_m));
if abs(alpha_m_mean - mean(alpha_m))<=0.01 && ...
   abs(pc_m_mean - mean(pc_m))<=0.01;
        bool = 1;
        sol_pc = mean(pc_m);
        sol_alpha = mean(alpha_m);
else
    alpha_m_mean = mean(alpha_m);
    pc_m_mean = mean(pc_m);
    alpha = linspace(0.9*min(alpha_m),1.1*max(alpha_m),50);
    pc = linspace(0.9*min(pc_m),1.1*max(pc_m),50);
end
end
% contourf(pc,rc,area);
% set(gca,'XTick',0:pi/10:2*pi)
% set(gca,'XTickLabel',{'0','pi/10','2*pi/10','3*pi/10','4*pi/10','5*pi/10','6*pi/10'...
%                        '7*pi/10','8*pi/10','9*pi/10','10*pi/10','11*pi/10','12*pi/10',...
%                        '13*pi/10','14*pi/10','15*pi/10','16*pi/10','17*pi/10','18*pi/10','19*pi/10'})