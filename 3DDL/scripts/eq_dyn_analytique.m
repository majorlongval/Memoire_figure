%% Script pour déterminer les équations de la dynamique du système.
clear all; close all; clc;

% Déclaration des variables
syms r R alpha x y z r_c phi_c h_c ddx ddy ddz g m tx ty tz Mx My Mz c_x c_y c_z Tz real 

% Création des vecteurs
p = [x;y;z];

Qs = [cos(2*pi/3) -sin(2*pi/3) 0;
    sin(2*pi/3)  cos(2*pi/3) 0;
    0                     0  1];

a11 = r*[0;1;0];
a12 = -a11;

a21 = Qs*a11;
a22 = -a21;

a31 = Qs*a21;
a32 = -a31;

b1m = R*[cos(alpha);sin(alpha);0];
b2m = Qs*b1m;
b3m = Qs*b2m;

c = [r_c*cos(phi_c);r_c*sin(phi_c);h_c];


e1 = p-b1m;
e2 = p-b2m;
e3 = p-b3m;

delta11 = cross(a11-c,e1);
delta12 = cross(a12-c,e1);
delta21 = cross(a21-c,e2);
delta22 = cross(a22-c,e2);
delta31 = cross(a31-c,e3);
delta32 = cross(a32-c,e3);

M = [e1 e1 e2 e2 e3 e3;
    delta11 delta12 delta21 delta22 delta31 delta32];
Msing = simplify(subs(M,alpha,pi/2));

detM = simplify(det(M'));

adjM = adjoint(M);

gamma = [-m*(ddx)+tx;-m*ddy+ty;m*(g-ddz)+tz;+Mx;+My;+Mz];

cond = simplify(adjM/(9*R^2*r^2*z^2));








M_aug = [M,[tx;ty;tz;Mx;My;Mz]];

coeffX = sym(zeros(6,6));
coeffY = sym(zeros(6,6));
coeffZ = sym(zeros(6,6));
coeffK = sym(zeros(6,6));

for i=1:6
    for j=1:6
        [T,C] = coeffs(cond(i,j),[x,y,z]);
        for k =1:length(C)
           if boolean(C(k) == x)
               coeffX(i,j) = T(k);
           elseif boolean(C(k) == y)
               coeffY(i,j) = T(k);
           elseif boolean(C(k) == z) 
               coeffZ(i,j) = T(k);
           elseif boolean(C(k) == 1)
               coeffK(i,j) = T(k);
           end
        end
    end
end
% 
% for i =1:6
%     for j =1:6
%         CoeffXlatex{i,j} = latex(coeffX(i,j));
%         CoeffYlatex{i,j} = latex(coeffY(i,j));
%         CoeffZlatex{i,j} = latex(coeffZ(i,j));
%         CoeffKlatex{i,j} = latex(coeffK(i,j));
%     end
% end
% 
% 
% filecoeff_X = fopen('coeff_X.txt','w');
% [rows,cols]=size(CoeffXlatex);
%  for r=1:rows
%      for c=1:cols
%         fprintf(filecoeff_X,'%s,    ',CoeffXlatex{r,c});
%      end
%      fprintf(filecoeff_X,'\n');
%  end
%  
% 
% filecoeff_Y = fopen('coeff_Y.txt','w');
% [rows,cols]=size(CoeffYlatex);
%  for r=1:rows
%      for c=1:cols
%         fprintf(filecoeff_Y,'%s,    ',CoeffYlatex{r,c});
%      end
%      fprintf(filecoeff_Y,'\n');
%  end
%  
% 
% filecoeff_Z = fopen('coeff_Z.txt','w');
% [rows,cols]=size(CoeffZlatex);
%  for r=1:rows
%      for c=1:cols
%         fprintf(filecoeff_Z,'%s,    ',CoeffZlatex{r,c});
%      end
%      fprintf(filecoeff_Z,'\n');
%  end
%  
%  
% filecoeff_K = fopen('coeff_K.txt','w');
% [rows,cols]=size(CoeffKlatex);
%  for r=1:rows
%      for c=1:cols
%         fprintf(filecoeff_K,'%s,    ',CoeffKlatex{r,c});
%      end
%      fprintf(filecoeff_K,'\n');
%  end
%  
%  
%  
% for i =1:6
%     Mt(i,:) = [coeffX(i,2)*ty+coeffX(i,3)*Tz+coeffX(i,4)*Mx,...
%                coeffY(i,1)*tx+coeffY(i,3)*Tz+coeffX(i,5)*My,...
%                coeffZ(i,1)*tx+coeffZ(i,2)*ty+coeffZ(i,6)*Mz];
% end

% Mt = [e1 e2 e3];
% detMt = det(Mt);
%
%
% cond1 = simplify(cross(a11+a21,b1m+b2m));
% cond2 = simplify(cross(a11+a22,b1m+b2m));
% cond3 = simplify(cross(a11+a31,b1m+b3m));
% cond4 = simplify(cross(a11+a32,b1m+b3m));
%
% cond5 = simplify(cross(a12+a21,b1m+b2m));
% cond6 = simplify(cross(a12+a22,b1m+b2m));
% cond7 = simplify(cross(a12+a31,b1m+b3m));
% cond8 = simplify(cross(a12+a32,b1m+b3m));
%
%
% zeta11x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,0));
% zeta11y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,0));
%
% zeta12x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,pi/2));
% zeta12y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,pi/2));
%
% zeta13x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,pi));
% zeta13y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,pi));
%
% zeta14x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,3*pi/2));
% zeta14y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
% zeta21x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,0));
% zeta21y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,0));
%
% zeta22x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,pi/2));
% zeta22y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,pi/2));
%
% zeta23x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,pi));
% zeta23y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,pi));
%
% zeta24x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a11(1))/(b2m(1)-b1m(1)),alpha,3*pi/2));
% zeta24y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a11(2))/(b2m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
%
% zeta31x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,0));
% zeta31y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,0));
%
% zeta32x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,pi/2));
% zeta32y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,pi/2));
%
% zeta33x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,pi));
% zeta33y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,pi));
%
% zeta34x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,3*pi/2));
% zeta34y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
%
% zeta41x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,0));
% zeta41y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,0));
%
% zeta42x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,pi/2));
% zeta42y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,pi/2));
%
% zeta43x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,pi));
% zeta43y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,pi));
%
% zeta44x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a11(1))/(b3m(1)-b1m(1)),alpha,3*pi/2));
% zeta44y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a11(2))/(b3m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
%
% zeta51x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,0));
% zeta51y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,0));
%
% zeta52x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,pi/2));
% zeta52y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,pi/2));
%
% zeta53x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,pi));
% zeta53y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,pi));
%
% zeta54x = simplify(subs((b2m(1)+a21(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,3*pi/2));
% zeta54y = simplify(subs((b2m(2)+a21(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
% zeta61x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,0));
% zeta61y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,0));
%
% zeta62x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,pi/2));
% zeta62y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,pi/2));
%
% zeta63x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,pi));
% zeta63y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,pi));
%
% zeta64x = simplify(subs((b2m(1)+a22(1)-b1m(1)-a12(1))/(b2m(1)-b1m(1)),alpha,3*pi/2));
% zeta64y = simplify(subs((b2m(2)+a22(2)-b1m(2)-a12(2))/(b2m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
% zeta71x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,0));
% zeta71y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,0));
%
% zeta72x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,pi/2));
% zeta72y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,pi/2));
%
% zeta73x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,pi));
% zeta73y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,pi));
%
% zeta74x = simplify(subs((b3m(1)+a31(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,3*pi/2));
% zeta74y = simplify(subs((b3m(2)+a31(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,3*pi/2));
%
%
%
%
%
% zeta81x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,0));
% zeta81y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,0));
%
% zeta82x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,pi/2));
% zeta82y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,pi/2));
%
% zeta83x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,pi));
% zeta83y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,pi));
%
% zeta84x = simplify(subs((b3m(1)+a32(1)-b1m(1)-a12(1))/(b3m(1)-b1m(1)),alpha,3*pi/2));
% zeta84y = simplify(subs((b3m(2)+a32(2)-b1m(2)-a12(2))/(b3m(2)-b1m(2)),alpha,3*pi/2));
