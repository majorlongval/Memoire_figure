%% Script pour déterminer les équations de la dynamique du système.
clear all; close all; clc;

% Déclaration des variables
syms r R alpha x y z r_c phi_c h_c ddx ddy ddz g m tx ty tz Mx My Mz c_x c_y c_z Tz real 

% Création des vecteurs
p = [x;y;z];

Qs = [cos(2*pi/3) -sin(2*pi/3) 0;
    sin(2*pi/3)  cos(2*pi/3) 0;
    0                     0  1];
a11 = r*[1/2;sqrt(3)/2;0];
a12 = r*[-1/2;sqrt(3)/2;0];

a21 = Qs*a11;
a22 = Qs*a12;

a31 = Qs*a21;
a32 = Qs*a22;

R1  = R*[0;1;0];
R2  = Qs*R1;
R3  = Qs*R2;


