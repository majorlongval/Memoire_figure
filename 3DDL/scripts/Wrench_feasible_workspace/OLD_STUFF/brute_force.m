%% Brute force method to calculate the wrench feasible workspace
clear all;  clc;

addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
% Setting the arbitrary values
R     = 1;      % m
r     = 0.2;    % m
alpha = 0;   % rad
rc    = 0;   % m
phic  = pi/6;   % rad
hc    = 0;    % m


% Setting the maximum z
zmax = 10;


% Maximal and minimal cable tension
fmax = 50; % N
fmin = 0.2;   % N

% Range of values for the external wrench and the mass
Mx = 0;%[-0.4, -0.2];%[-0.02 0.02]; % Nm
My = 0;%[-0.01 0.05]; % Nm
Mz = 0;%[-0.1 0.1]; % Nm
fx = [-1,1]; % N
fy = [-1,1]; % N
fz = 0;%[-10  10]; % N
m  = 3; %kg

% Span of the z axis incrementation
z_span = linspace(0,zmax,30);


% Discretising the values to check
x = linspace(-1,1,40);
y = linspace(-1,1,40);
pos = [];
for i =1:length(z_span)
    for j = 1:length(y)
        for k = 1:length(x)
            pos = [pos,[x(k);y(j);z_span(i)]];
        end
    end
end



% Determining which points are part of the WFW
good_points = [];
for i =1:length(pos(1,:))
    return_bool = check_cond(pos(:,i),R,r,alpha,rc,phic,hc,fx,fy,fz,Mx,My,Mz,m,fmin,fmax);
    if return_bool == 1
        good_points = [good_points,pos(:,i)];
    end
end



plot3(good_points(1,:),good_points(2,:),good_points(3,:),'*g');
set(gca,'Zdir','reverse');
set(gca,'Xdir','reverse');
grid on;
xlabel('XXX');
ylabel('YYY');
zlabel('ZZZ');
title({'TITLE1','TITLE2'});
view([-79,34]);









%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function return_bool = check_cond(pos,R,r,alpha,rc,phic,hc,fx,fy,fz,Mx,My,Mz,m,fmin,fmax)
wrench_comb = [];
for a = 1:length(fx)
    for b =1:length(fy)
        for c = 1:length(fz)
            for d = 1:length(Mx)
                for e = 1:length(My)
                    for f = 1:length(Mz)
                        wrench_comb = [wrench_comb;...
                            fx(a),fy(b),fz(c),Mx(d),My(e),Mz(f)];
                    end
                end
            end
        end
    end
end
for k =1:length(wrench_comb(:,1))
tx = wrench_comb(k,1);
ty = wrench_comb(k,2);
tz = wrench_comb(k,3);
Mx = wrench_comb(k,4);
My = wrench_comb(k,5);
Mz = wrench_comb(k,6);
g =9.81;
Tz = tz + m*g;

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, rc, phic, hc);

rho = calc_rho(pos, R, alpha);
roundn = @(x,n) round(x.*10.^n)./10.^n;
for i = 1:6
    A(i) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
    B(i) = tx*ay(i) + Tz*cy(i) + My*ey(i);
    C(i) = tx*(az(i)*pos(3)+ak(i))+ty*(bz(i)*pos(3)+bk(i))+Tz*ck(i)+...
        Mx*dk(i)+My*ek(i)+Mz*(hz(i)*pos(3)+hk(i))-(6*R*r*pos(3)*cos(alpha)*fmin/rho(i));
    if (roundn(A(i)*pos(1)+B(i)*pos(2)+C(i),4)>=0)
        list_elem(i) = 1;
    else
        list_elem(i) = 0;
    end
end
for i = 7:12
    A(i,1) = -A(i-6);
    B(i,1) = -B(i-6);
    C(i,1) = -(tx*(az(i-6)*pos(3)+ak(i-6))+ty*(bz(i-6)*pos(3)+bk(i-6))+Tz*ck(i-6)+...
        Mx*dk(i-6)+My*ek(i-6)+Mz*(hz(i-6)*pos(3)+hk(i-6)))+(6*R*r*pos(3)*cos(alpha)*fmax/rho(i-6));
    if (roundn(A(i)*pos(1)+B(i)*pos(2)+C(i),4)>=0)
        list_elem(i) = 1;
    else
        list_elem(i) = 0;
    end
end

if all(list_elem(:) == 1)
    rb(k) = 1;
else
    rb(k) = 0;
end
end

if all(rb(:) ==1)
    return_bool = 1;
else
    return_bool = 0;
end

end


function rho = calc_rho(pos, R, alpha)
Qs = [cos(2*pi/3) -sin(2*pi/3) 0;...
    sin(2*pi/3) cos(2*pi/3) 0;
    0 0 1];

RR(:,1) = R*[cos(alpha);sin(alpha);0];
RR(:,2) = RR(:,1);
RR(:,3) = Qs*RR(:,2);
RR(:,4) = Qs*RR(:,2);
RR(:,5) = Qs*RR(:,4);
RR(:,6) = Qs*RR(:,4);

for i =1:6
    rho(i) = norm(pos-RR(:,i));
end
end


