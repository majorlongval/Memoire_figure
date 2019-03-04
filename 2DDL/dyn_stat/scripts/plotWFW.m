%% Script to calculate the WFW for the planar mechanism
clear all; close all; clc;
addpath('functions');


x = 5; % m
y = 0; % m
% Defining a single geometry
geo = pack_geo(0,0,0,2,0.5);

% Defining the  two sets of wrenchs
rangefx1 = [-0.5,-0.1];
rangefy1 = [-0.1,0.0];
rangeMe1 = [-0.1,0.1];

rangefx2 = [-0.5,-0.1];
rangefy2 = [-0.5,0.5];
rangeMe2 = [-0.1,0.1];

[var11,var12,var13,var14] = Effort_poly(x, y, geo,rangefx1,rangefy1,rangeMe1);

[var21,var22,var23,var24] = Effort_poly(x, y, geo,rangefx2,rangefy2,rangeMe2);


saveas(var14,'good_WFW.svg');
saveas(var24,'bad_WFW.svg');