%% Script pour calculer la position des oeillets. 
clear all; close all; clc;

L = 2110; % mm longueur d'un côté du triangle
r_ext = sqrt(3)*L/6; % mm rayon du cercle inscrit

alpha = pi/6;

d_M = tan(alpha)*r_ext;

pos_mot_2= [r_ext;d_M;0];

attache_21 =  [0;0.2;0];
attache_22 =  [0;-0.2;0];


cs = cos(2*pi/3);
ss = sin(2*pi/3);

Qs = [cs -ss 0;
      ss  cs 0;
       0   0 1];
pos_mot_3 = Qs*pos_mot_2;
attache_31 = Qs*attache_21;
attahce_32 = Qs*attache_22;

pos_mot_1 = Qs*pos_mot_3;
attache_11 = Qs*attache_31;
attache_12 = Qs*attahce_32;

ini_hauteur = 3204.4; % mm hateur initiale du mécanisme

pos_ini = [0; 0; 3204.4];

cable_2_l_ini = norm(pos_ini-pos_mot_2);
cable_3_l_ini = norm(pos_ini-pos_mot_3);
cable_1_l_ini = norm(pos_ini-pos_mot_1);
