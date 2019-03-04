%% Optimisation des paramètres alpha et pc
clear all;close all;clc;

% Adding the path to the functions 
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');


% The size of the big and small circles
R = 1;
r = 0.2;

% The span of rc on which I want to optimise
rc = linspace(0,r,100);


% The bounds of the alpha and phic parameters
lb = [0 0 ];
ub = [pi/2, (pi/3)];

% Vectors containing initial guesses for the location of the maximums.

alpha_0 = pi/6;
pc_0 = pi/4;
% Calcul des valeurs qui min 
x0 = [alpha_0,pc_0];
[pc_sol, alpha_sol] = opti_max_dETS_phic_alpha(rc, R, r, ub, lb, x0);




% figure
% for i =1:10
%     plot(rc'.*cos(pc_sol_vect(:,i)),rc'.*sin(pc_sol_vect(:,i)),'.r');
%     hold on;
% end


function [pc_sol, alpha_sol] = opti_max_dETS_phic_alpha(rcspan, R, r,  ub, lb, x0)
        fun = @(x)mean_calc_dETS_SW_new(x(1),R,r,rcspan,x(2));
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        nonlcon = [];
        %options = optimoptions('fmincon','Algorithm','sqp');
        x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
        alpha_sol = x(1);
        pc_sol = x(2);
end

function mean_dETS = mean_calc_dETS_SW_new(alpha,R,r,rcspan,phic)
    for i =1:length(rcspan)
        dETS(i) = calc_dETS_SW_new(alpha,R,r,rcspan(i),phic);
    end
    mean_dETS = mean(dETS);
end