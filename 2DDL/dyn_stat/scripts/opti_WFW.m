% Script to find the optimal value of S_WFW
clear all;close all; clc;
L = 5;
l = 1;
m = 1;
Tx = [-0.04,0.06]*9.81*m;%[-0.1 0.1];
Ty = [-0.05,0.05]*9.81*m;
Me = [-0.05, 0]*9.81*m;



Stest = Calc_S_WFW(Tx,Ty,Me,m,L,l,0,0,0);



fun = @(x)-Calc_S_WFW(Tx,Ty,Me,m,L,l,x(1),x(2),x(3));

x0 = [0,0,0];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-l,-l,-l];
ub = [l,l,l];
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
