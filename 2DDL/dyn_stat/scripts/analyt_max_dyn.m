%% test max @ w= sqrt(T/mxc)
clear all; close all;clc;
syms L l ay cy cx x y g m A B C D tx ty Me alpha  r1 r2  y0 y1 x0  dx x1 dy real 

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

Psi1 = A1+B1*y0;
Psi2 = A2+B2*y0;
Psi3 = A3+B3*y0;

Omega1 = C1+B1*x0;
Omega2 = C2+B2*x0;
Omega3 = C3+B3*x0;


% eq = simplify(((Psi3*tx-B3*ty*x)*(Psi2*tx-B2*ty*x+D2*Me))-((Psi2*tx-B2*ty*x)*(Psi3*tx-B3*ty*x+D3*Me)));
% 
% 
% tysub = tx*(L+y+ay)/x;
% 
% eq2 = simplify((D3*Me)/(Psi3*tx-B3*tysub*x));


%dx = x1-x0;
%dy = y1-y0;

xxs = r2*dx;
yys = r2*dy;
xs = x0+r1*dx;
ys = y0+r1*dy;

cond = (m*(g-xxs)+tx)*(A+B*ys)+(m*yys-ty)*(C+B*xs)+D*Me;


cond 