%% Script to analyse the influence of the parameters on the area of SW
clear all;close all;clc;
alpha =pi/8;
R = 1;
r = 0.2;
rc = linspace(0,r,100);
pc = linspace(0,2*pi,100);

for i =1:length(rc)
    for j =1:length(pc)
        area(i,j) = calc_area_SW_new(alpha,R,r,rc(i),pc(j));
    end
end

contourf(pc,rc,area);
set(gca,'XTick',0:pi/10:2*pi) 
set(gca,'XTickLabel',{'0','pi/10','2*pi/10','3*pi/10','4*pi/10','5*pi/10','6*pi/10'...
                       '7*pi/10','8*pi/10','9*pi/10','10*pi/10','11*pi/10','12*pi/10',...
                       '13*pi/10','14*pi/10','15*pi/10','16*pi/10','17*pi/10','18*pi/10','19*pi/10'}) 