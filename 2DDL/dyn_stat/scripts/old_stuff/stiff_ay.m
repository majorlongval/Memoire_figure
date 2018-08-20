%% Script to determine best ay parameters for average cx and cy
clear all; close all; clc;
addpath('functions');
L = 5;
l = 2;
ay = linspace(-l,l,51);
cx = linspace(-l,l,51);
cy = linspace(-l,l,51);


% The height at which the test is done :
x = L;

RLv = zeros(1,51);
RL  = zeros(51,51,11);
for i = 1:length(ay)
    for j = 1:length(cy)
        for k = 1:length(cx)
            geo = pack_geo(ay(i),cy(j),cx(k),L,l);
            ymin = y_min_stat(geo);
            yv = linspace(ymin,L,50);
            for m =1:length(yv)
                RLv(m) = rigid_lat([x,yv(m)],geo);
            end
            RL(j,k,i) = mean(RLv);
        end
    end
    RL_ay = RL(:,:,i);
    RL_ay_mean(i) = mean(RL_ay(:));
end

plot(ay,RL_ay_mean);