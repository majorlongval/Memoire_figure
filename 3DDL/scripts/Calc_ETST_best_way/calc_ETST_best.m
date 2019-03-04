clear all; close all; clc;
%% script to calculate the volume of the ETST
g = 9.81;
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions');
% Setting the arbitrary values
R     = 1;      % m
r     = 0.2;    % m
alpha = 0;      % rad
rc    = 0.1;      % m
phic  = pi/6;   % rad
hc    = 0;      % m
m     = 1;      % kg
Tx    = 1;      % N
Ty    = -1;      % N
tz    = 0;      % N
Tz    = m*g+tz; % N
Mx    = 0.1;      % Nm
My    = 0;      % Nm
Mz    = -0.1;      % Nm

% Calculating the coefficients
[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, rc, phic, hc);

% Setting the maximum z
zmax = 10;
zmin = 1;
% Setting the plane equations

for i =1:6
    Ap(i) = Ty*bx(i)+Tz*cx(i)+Mx*dx(i);
    Bp(i) = Tx*ay(i) + Tz*cy(i) + My*ey(i);
    Cp(i) = Tx*az(i) + Ty*bz(i) + Mz*hz(i);
    Dp(i) = Tx*ak(i) + Ty*bk(i) + Tz*ck(i)+Mx*dk(i)+My*ek(i)+Mz*hk(i);
end
Ap(7) = 0; Bp(7) = 0; Cp(7) = 1; Dp(7) = zmin;
Ap(8) = 0; Bp(8) = 0; Cp(8) = -1; Dp(8) = zmax;

Mp = [Ap',Bp',Cp',Dp'];

% Doing all the combinations and the solutions 

v = 1:size(Mp,1);
a = nchoosek(v,3);
B = zeros(4,3,length(a));
intersection = [];
a_inter  = [];
warning('off','MATLAB:singularMatrix');
for i = 1:length(a)
    B(:,:,i) = Mp(a(i,:),:)';
    C(:,i) = B(1:3,:,i)'\-B(4,:,i)';
    bool = true;
    for j =1:8
        if ~ismember(j,a(i,:))
            temp = Mp(j,:)*[C(:,i);1] > 0;
            bool = bool*temp;
        end
    end
    if bool == true
        intersection = [intersection, C(:,i)];
        a_inter = [a_inter;a(i,:)];
    end
end
X = intersection';
DTo = delaunayTriangulation(X);
[K,v] = convexHull(DTo);
% plot3(intersection(1,:),intersection(2,:),intersection(3,:),'*r');
tri_surf_h = trisurf(K,X(:,1),X(:,2),X(:,3));
tri_surf_h.FaceColor = 'r';
tri_surf_h.FaceAlpha = 0.5;









% %Calculating the  base of ETST at z = 0;
% for i =1:6
%     A0(i) = Ty*bx(i) + Tz*cx(i) + Mx*dx(i);
%     B0(i) = Tx*ay(i) + Tz*cy(i) + My*ey(i);
%     C0(i) = Tx*ak(i)+Ty*bk(i)+Tz*ck(i)+...
%         Mx*dk(i)+My*ek(i)+Mz*hk(i);
%     Line(i,:) = [A0(i) B0(i) C0(i)];
%     Ap(i) = Ty*bx(i)+Tz*cx(i)+Mx*dx(i);
%     Bp(i) = Tx*ay(i) + Tz*cy(i) + My*ey(i);
%     Cp(i) = Tx*az(i) + Ty*bz(i) + Mz*hz(i);
%     Dp(i) = Tx*ak(i) + Ty*bk(i) + Tz*ck(i)+Mx*dk(i)+My*ek(i)+Mz*hk(i);
% end
% 
% intersection = cell(6,2);
%  for i =1:5
%      for j =i+1:6
%          p_inter = [A0(i) B0(i);A0(j) B0(j)]\-[C0(i);C0(j)];
%          temp_i = intersection{i,2};
%          temp_i = [temp_i,p_inter];
%          intersection{i,2} = temp_i;
%          temp_i2 = intersection{i,1};
%          temp_i2 = [temp_i2,j];
%          intersection{i,1} =temp_i2;
%          temp_j = intersection{j,2};
%          temp_j = [temp_j,p_inter];
%          intersection{j,2} = temp_j;
%          temp_j2 = intersection{j,1};
%          temp_j2 = [temp_j2,i];
%          intersection{j,1} = temp_j2;
%      end
%  end
% vert = [];
% vert_l = [];
% 
%  for i =1:6
%     temp = intersection{i,2};
%     [Xp, Yp] = vec_acw_order_lin(temp(1,:)',temp(2,:)');
%     
%     for j = 1:length(Xp)-1
%         mid = 0.5*[Xp(j)+Xp(j+1),Yp(j)+Yp(j+1)];
%     bool =1;
%     for k =1:length(Line)
%         if k ~= i
%             tbool = A0(k)*mid(1)+B0(k)*mid(2)+C0(k)>0;
%             bool = bool*tbool;
%         end
%     end
%     if bool == 1
%         vert = [vert;Xp(j),Yp(j);Xp(j+1),Yp(j+1)];
%         tempb = intersection{i,1};
%         vert_l = [vert_l;i,tempb(j);i,tempb(j+1)];
%     end
%     end
%  end
% sol = {};
%  for i =1:length(vert_l(:,1))
%      i1 = vert_l(i,1);
%      i2 = vert_l(i,2);
%      for j =1:6
%          if ~any(vert_l(i,:)==i)
%             temp = [Ap(i1) Bp(i1) Cp(i1);...
%                     Ap(i2) Bp(i2) Cp(i2);...
%                     Ap(j)  Bp(j)  Cp(j) ]\[Dp(i1);Dp(i2);Dp(j)];
%             sol = {sol,temp};
%          end
%      end
%  end
