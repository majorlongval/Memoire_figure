function [volume,K,corners] = WFW_Calculator(R,r,alpha,rc,pc,hc,m,wc,zr)
%WFW_Calculator Calculates the corners of the WFW between zmin and zmax for
% a given geometrical arrangement. It also gives the volume on the given
% height
%   coners gives a struct of which has the following structure:
% Plane 1, Plane 2, Plane 3, x, y, z
%   Volume is the volume of the WFW on the given height




% Calculating all the possible planes
planes = [];
for k =1:length(wc(:,1))
    planes = [planes;Calc_planes(alpha,R,r,rc,pc,hc,m,...
        wc(k,:))];
end
planes = [planes;0 0 -1,zr(2);0 0 1 -zr(1)];

% Defining the list of all projection lines at the plane z =zmin
L_l_i = [];
for i =1:length(wc(:,1))
    L_l_i = [L_l_i;Calc_lines(alpha, R, r, rc, pc, ...
        hc,m,zr(1),wc(i,:))];
end

% Determining the points on the z=zmin plane which are good.

pii = [];
for i =1:length(L_l_i(:,1))-1
    for j=i+1:length(L_l_i(:,1))
        temp = Calc_intersection(L_l_i(i,:), L_l_i(j,:));
        bool = 1;
        for k =1:length(L_l_i(:,1))
            if k ~= i && k~= j
                verif = L_l_i(k,1)*temp(1)+L_l_i(k,2)*temp(2)+L_l_i(k,3) >= 0;
                bool = bool*verif;
            end
        end
        if bool ==1
            pii = [pii;i,j,194,temp,zr(1)];
            
        end
    end
end

% Determining all_the points on the line and determining the shortest
% segment
for i =1:length(pii(:,1))
    dist = [];
    for j =1:length(planes(:,1))
        if ~ismember(j,pii(i,1:3))
            x = calc_plane_intersection(planes(pii(i,1),:),...
                planes(pii(i,2),:),...
                planes(j,:));
            if ~isempty(x) && x(3)> pii(i,6)
                dist = [dist;j,dist_2_points(pii(i,4:6),x)];
            end
        end
    end
    dist = sortrows(dist,2,'ascend');
    short = dist(1,:);
    pii = [pii; pii(i,1:2),short(1),calc_plane_intersection(planes(pii(i,1),:),...
                                           planes(pii(i,2),:),...
                                           planes(short(1),:))'];
end

corners = pii;

[K,volume] = convhulln(corners(:,4:6));





%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Planes = Calc_planes(alpha, R, r, r_c, phi_c, ...
                                         h_c,m,wrench_comb)
%Calc_planes calculates the equations of planes as a function of the geometric parameters
%in the form Ax+By+Cz+D = 0
g =9.81;
tx = wrench_comb(1);
ty = wrench_comb(2);
tz = wrench_comb(3);
Mx = wrench_comb(4);
My = wrench_comb(5);
Mz = wrench_comb(6);

Tz = tz + m*g;

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, r_c, phi_c, h_c);
for i = 1:6
    A(i,1) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
    B(i,1) = tx*ay(i) + Tz*cy(i) + My*ey(i);
    C(i,1) = tx*az(i)+ty*bz(i)+Mz*hz(i);
    D(i,1) = tx*ak(i)+ty*bk(i)+Tz*ck(i)+Mx*dk(i)+My*ek(i)+Mz*hk(i);

end
Planes = [A, B, C, D];
end

function X = calc_plane_intersection(p1,p2,p3)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
M =[p1(1),p1(2),p1(3);...
    p2(1),p2(2),p2(3);...
    p3(1),p3(2),p3(3)];
if rank(M) == 3
    d = -[p1(4);p2(4);p3(4)];
    X = M\d;
else
    X = [];
    
    
    
end


end

function dist = dist_2_points(p1,p2)
    dist = sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2+(p1(3)-p2(3))^2);
end

function List_lines = Calc_lines(alpha, R, r, r_c, phi_c, ...
                                         h_c,m,z0,wrench_comb)
%CALC_LINES Calculates the A B and C coefficients of the lines associated
%with a set of wrench 
%   List_lines is a matrix of the form [A B C] where A B C are column
%   vectors of length 6.
g =9.81;
tx = wrench_comb(1);
ty = wrench_comb(2);
tz = wrench_comb(3);
Mx = wrench_comb(4);
My = wrench_comb(5);
Mz = wrench_comb(6);

Tz = tz + m*g;

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, r_c, phi_c, h_c);
for i = 1:6
    A(i,1) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
    B(i,1) = tx*ay(i) + Tz*cy(i) + My*ey(i);
    C(i,1) = tx*(az(i)*z0+ak(i))+ty*(bz(i)*z0+bk(i))+Tz*ck(i)+...
        Mx*dk(i)+My*ek(i)+Mz*(hz(i)*z0+hk(i));
end
List_lines = [A,B,C];
end

function p_inter = Calc_intersection(D1, D2)
%CALC_INTERSECTION calculates all the intersections between D1 and D2
%   D1 and D2 and lines that are given by Di = [Ai, Bi, Ci] such that
%   Aix+Biy+Ci = 0;
p_inter = (1/(D1(2)*D2(1)-D1(1)*D2(2)))*...
          [D1(3)*D2(2)-D1(2)*D2(3), D1(1)*D2(3)-D2(1)*D1(3)];
end

end
