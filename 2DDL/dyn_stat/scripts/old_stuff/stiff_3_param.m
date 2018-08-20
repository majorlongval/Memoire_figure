%% Script to determine the latteral stifness as a function of all three param.
clear all; close all; clc;
addpath('functions');
L = 5;
l = 1;
ay = linspace(-l,l,51);
cy = ay;
cx = linspace(-l,l,51);


% The height at which the test is done :
x = L;
RLv = zeros(1,51);
RL  = zeros(51,51,11);
for i = 1:length(cx)
    for j = 1:length(ay)
        for k = 1:length(cy)
            geo = pack_geo(ay(j),cy(k),cx(i),L,l);
            ymin = y_min_stat(geo);
            yv = linspace(ymin,L,50);
            for m =1:length(yv)
                RLv(m) = rigid_lat([x,yv(m)],geo);
            end
            RL(j,k,i) = mean(RLv);
        end
    end
    RLtest = RL(:,:,i);
   [RLmax(i),idx(i)] = max(RLtest(:));
end