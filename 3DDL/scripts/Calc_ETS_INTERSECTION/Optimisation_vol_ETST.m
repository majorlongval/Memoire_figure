%% Optimisation des paramètres géométriques pour maximiser le volume de l'ETST
clear all; clc;

% Adding the path to the functions 
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

% The size of the big and small circles
R = 1;    %m
r = 0.2;  %m
m = 1;
% The bounds of the parameters
lb = [ 0 0.0001 ];
ub = [ pi/3 pi/2];


% Vector of initial guess

x0 = (ub + lb)/2;


rc = 0.02;
hc = 0.1;


% Creating the possible values of the forces and moments
fx = linspace(-1,1,5);
fy = linspace(-1,1,5);
fz = 0.3;
Mx = 0.1;
My = 0.1;
Mz = 0.1;
comb_force = [];
for a = 1:length(fx)
    for b = 1:length(fy)
        for c = 1:length(fz)
            for d = 1:length(Mx)
                for e = 1:length(My)
                    for f = 1:length(Mz)
                        comb_force = [comb_force;...
                            fx(a), fy(b), fz(c), Mx(d), My(e), Mz(f)];
                    end
                end
            end
        end
    end
end

x_sol = opti_max_surf_phic_alpha(R, r,  ub, lb, x0,rc,hc,comb_force,m);
 
function x_sol = opti_max_surf_phic_alpha(R, r,  ub, lb, x0,rc,hc,comb_force,m)
        fun = @(x)-mean_calc_vol_SW_Inter(x,rc,hc,R,r,comb_force,m);
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        nonlcon = [];
        %options = optimoptions('fmincon','Algorithm','sqp');
        options = optimoptions('fmincon','Display','iter');
        x_sol = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end



function mean_vol = mean_calc_vol_SW_Inter(x,rc,hc,R,r,comb_force,m);
    vol_tot = 0;
    for i =1:length(comb_force(:,1))
        vol_tot = vol_tot + calc_vol_INTER(x(2), R, r, rc, x(1), ...
    hc, comb_force(i,1), comb_force(i,2), comb_force(i,3), m, comb_force(i,4),...
    comb_force(i,5), comb_force(i,6),10,50);
%     display(i);
    end
    mean_vol = vol_tot/i;
%     display('mean_calculated');
end