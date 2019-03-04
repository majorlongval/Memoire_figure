%% Script for analytical formula of stifness
clear all; close all; clc;
syms L l cy ay cx x y g r1 r2  fx fy Me real

a1 = [0;ay];
a2 = [0;-l];
a3 = [0;l];

b1 = [0;-L];
b2 = [0;L-l];
b3 = [0;L+l];

pos = [x;y];

c = [cx;cy];

E = [0 -1;1 0];

gam = [g;0;0];

e1 = (pos+a1-b1)/r1;
e2 = (pos+a2-b2)/r2;
delta1    = (a1-c)'*E*e1;
delta2    = (a2-c)'*E*e2;
delta3    = (a3-c)'*E*e2;



gamma = [g+fx; fy; Me];


detM = (a2-a3)'*E*e2*(e1'*E*e2)

M = [e1' delta1;
     e2' delta2;
     e2' delta3];
detMspec = simplify(det(subs(M,ay,l)));

M
% 
% sol = adjoint([e1 e2 e2;delta1 delta2 delta3])*gamma;
% % % Latteral stifness
% ey = [g+fx;fy];
% 
% M1 = [ey e2 e2;Me delta2 delta3];
% M2 = [e1 ey e2;delta1 Me delta3];
% M3 = [e1 e2 ey;delta1 delta2 Me];
% 
% sol1t =  collect(simplify(det(M1))==0,[fx,fy,Me,g]);
% sol2t = collect(simplify(det(M2))==0,[fx,fy,Me,g]);
% sol3t = collect(simplify(det(M3))==0,[fx,fy,Me,g]);





