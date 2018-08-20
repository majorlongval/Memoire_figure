%% This script plots the stifness on the xy cy plane for given (ay,cy,cx).

clear all; close all; clc;
addpath('functions');

% Setting the geometric parameters
L = 5;
l = 1;
cx = -0.1;
cy = 0.1;
ay = -0.1;

gp = struct('L',L,'l',l,'ay',ay,'cy',cy,'cx',cx);

B2   = ay-l;
B3   = -(ay+l);
A3 = ay*(L-cy)-2*L*cy-l*(ay+L);
A2 = -ay*(L-cy)+2*L*cy-l*(ay+L);
ymin  = max(-A3/B3,-A2/B2);
% Creating an xy discrete field
xpos = linspace(1,5,100);
ypos = linspace(ymin, L, 100);

% Pre-allocating space 
Stiff(:,:,1) = zeros(length(xpos),length(ypos));
Stiff(:,:,2) = zeros(length(xpos),length(ypos));

% Stiff(:,:,1) is for latteral stiff and Stiff(:,:,2) is for rotational

% Calculation for every point in the field
for i =1:length(xpos)
    for j =1:length(ypos)
            [Stiff(i,j,1),Stiff(i,j,2)] = Rigid([xpos(i),ypos(j)],gp);
        
    end
end

figure
surf(ypos,xpos,Stiff(:,:,1));

figure
surf(ypos,xpos,Stiff(:,:,2));
