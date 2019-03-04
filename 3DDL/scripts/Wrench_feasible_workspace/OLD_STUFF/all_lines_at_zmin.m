clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');

% Setting the arbitrary values
R     = 0.75;                 % m
r     = 0.2;                  % m
alpha = pi/6;                 % rad
rc    = 0.01;              % m
phic  = 0;                    % rad
hc    = 0;                    % m


% Setting the minimum and maximum heigth
zmin = 2;                     % m
zmax = 10;                    % m

% Range of values for the external wrench and the mass
Mx = [-1 1];%[-1 1];                       % Nm
My = [-1 1];                       % Nm
Mz = [-1 1];                       % Nm
fx = [-1 1];                  % N
fy = [-1 1];                       % N
fz = 0;%[-1 1];                       % N
m  = 5;                     % kg

% Defining all combinations of wrenches ( can define multiple fo mult.
% plot)
wc = wrench_comb(fx,fy,fz,Mz,My,Mz);




% Defining the list of all projection lines at the plane z =zmin
L_l_i = [];
for i =1:length(wc(:,1))
    L_l_i = [L_l_i;Calc_lines(alpha, R, r, rc, phic, ...
        hc,m,zmin,wc(i,:))];
end
figure;
for i =1:length(L_l_i(:,1))
    plot_line(L_l_i(i,1),L_l_i(i,2),L_l_i(i,3),-5,5);
    hold on;
end
axis([-0.8 0.8 -0.8 0.8]);



% Determining all good points of intersection
pii = [];
counter  = 0;
for i =1:length(L_l_i(:,1))
    for j=i+1:length(L_l_i(:,1))
        temp = Calc_intersection(L_l_i(i,:), L_l_i(j,:));
        bool = 1;
        for k =1:length(L_l_i(:,1))
            if k ~= i && k~= j
                verif = L_l_i(k,1)*temp(1)+L_l_i(k,2)*temp(2)+L_l_i(k,3) >= 0;
                bool = bool*verif;
                counter = counter+1;
            end
        end
        if bool ==1
            pii = [pii;i,j,194,temp,zmin];
            plot(temp(1),temp(2),'.r','MarkerSize',20);
            
        end
    end
end
grid on;
saveas(gca,'all_lines_at_zmin','svg');


function wc = wrench_comb(fx,fy,fz,Mx,My,Mz)
wc = [];
for a = 1:length(fx)
    for b =1:length(fy)
        for c = 1:length(fz)
            for d = 1:length(Mx)
                for e = 1:length(My)
                    for f = 1:length(Mz)
                        wc = [wc;...
                            fx(a),fy(b),fz(c),Mx(d),My(e),Mz(f)];
                    end
                end
            end
        end
    end
end
end

function plot_line(a,b,c,xmin,xmax)
f =  @(x) (-a*x-c)/b;
plot([xmin,xmax],[f(xmin),f(xmax)],'-b');
end