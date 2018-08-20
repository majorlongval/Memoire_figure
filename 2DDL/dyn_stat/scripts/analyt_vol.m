% Calcul_analyt du volume
clear all; close all;clc;
syms L l ay cy cx x y g m alpha real

A1 = -L;
A2 = ay*(cy-L)+2*L*cy-l*(ay+L);
A3 = ay*(L-cy)-2*L*cy-l*(ay+L);

B1 = 1;
B2 = ay-l;
B3 = -(ay+l);

C1   = 0;
C2   = cx*(2*L+ay);
C3   = -C2;
D1   = 0;
D2   = 2*L+ay;
D3   = -D2;

u(:,1) = [0;[-C1-B1*x D1; -C2-B2*x D2]\...
       (-(alpha+g)*[A1+B1*y;A2+B2*y])]-[-g; 0; 0];
u(:,2) = [0;[-C1-B1*x D1; -C3-B3*x D3]\...
       (-(alpha+g)*[A1+B1*y;A3+B3*y])]-[-g; 0; 0];
u(:,3) = [0;[-C2-B2*x D2; -C3-B3*x D3]\...
       (-(alpha+g)*[A2+B2*y;A3+B3*y])]-[-g; 0; 0];

V = det(u);
   

% Calcul inutile
M  = [A1+B1*y -(C2+B2*x) D1;
      A2+B2*y -(C2+B2*x) D2;
      A3+B3*y -(C3+B3*x) D3];

b = -m*g*[A1+B1*y;
        A2+B2*y;
        A3+B3*y];
    
sol = M\b;
        
      



% 
% % Calcul des ranges de dfy et dMe
% %dfy
% %cas 1 : 2>3
% dfy1 = (-g*(A1+B1*y)/(C1+B1*x))+(g*(A2+B2*y)/(C2+B2*x));
% dfy1int = int(dfy1,y,-L,L);
% dfy1moy = simplify(dfy1int/(2*L));
% dfy1moymoy = int(dfy1moy,cx,-l,l)/(2*l);
% %dfy1moymoymoy = int(dfy1moymoy,cx,-l,l)/(2*l);
% 
% 
% %cas 1 : 2<3
% dfy2 = (-g*(A1+B1*y)/(C1+B1*x))+(g*(A3+B3*y)/(C3+B3*x));
% dfy2int = int(dfy2,y,-L,L);
% dfy2moy = simplify(dfy2int/(2*L));
% 
% %dMe
% dMe = g*(A2+A3+y*(B2+B3))/D2;
% dMemoy = simplify(int(dMe,y,-L,L));
% 
