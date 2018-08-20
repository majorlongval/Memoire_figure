%% This script calculates the evolution of the surface of the SW
%{
To do this, the integral of the static workspace will be calculated for
each possible combination of the values (ax,cy,ay).

%}
clear all; close all; clc;
addpath('functions');
%% 1 - Creating the geometric parameters

% The base parameters
L = 5;
l = 1;


%% 3 - Creating the XY field for testing the WS surface
xmax = L;
xmin = 0.01; % So not to have singularities in the calculations;
ymin = -(L+l);
ymax = L;

Surfmax = (ymax-ymin)*(xmax-xmin); % The rectangle giving the max integ.

% Creating a vector along x to get the ymin curve
x = xmin:0.1:xmax;

%% 2 - Creating vectors for ax,ay,cy
ax = -l:0.1:l;
ay = ax;
cy = ax;
cx = 0; %Changes nothing in the static ws.

% A struct containing the base parameters
for i =1:length(ax)
    for j =1:length(ay)
        for k =1:length(cy)
            gp(k,j,i) = struct('L',L,'l',l,'cy',cy(k),...
                'ay',ay(j),'cx',cx,'ax',ax(i));
            for m = 1:length(x)
                ystat(m,k,j,i) = dom_stat(x(m),gp(k,j,i));
                if ~isreal(ystat(m,k,j,i))
                    ystat(m,k,j,i) = NaN;
                end
            end
        end
    end
end

for i =1:length(ax)
    for j =1:length(ay)
        for k =1:length(cy)
            ystat_integ(k,j,i) = ystat_surf_integ(x,ystat(:,k,j,i));
            surf_integ_rel(k,j,i) = (Surfmax-ystat_integ(k,j,i))/(Surfmax);
        end
    end
end


sliceomatic(surf_integ_rel,ax,ay,cy);


figure
yt = linspace(-5,5,length(x));
for i =1:length(x)
    taut(:,i) = tension_stat([x(end),yt(i)],gp(2,15,13));
end

plot(yt,taut);
grid minor;
hold on;


figure
for i = 1:length(x)
    ytest_stat(i)= dom_stat(x(i),gp(3,15,13));
end
plot(x,ytest_stat)



plot([ystat(end,7,11,1),ystat(end,7,11,1)],[-5,10]);
