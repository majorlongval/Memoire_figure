%% exemple problem
clear all; close all; clc;
addpath('functions');
g =9.81;

%masses
mR = 1; %kg
m = 1; %kg
beta = mR/(m+mR);

% dim of R
b0 = 0.5; %m
b1 = 1;   %m
b2 = 1;   %m

% Position of C
cR = [0.8;-b0];
cE = beta*cR;

% Geometric parameters
L  = 5;   %m
l  = 1;   %m
ay = linspace(-l,l-0.01,100); %m

for i =1:length(ay)
    geo(i) = pack_geo(ay(i),cE(2),cE(1),L,l);
    [ay(i), cE(2), cE(1), L, l,A(:,i), B(:,i), C(:,i)] = unpack_geo(geo(i));
    D(:,i) = C(:,i)./cE(1);
end

% Force externe à appliquer

f = linspace(0,10,100); %N

% Torseurs total externe
for i =1:length(f)
    T(:,i) = [mR*g;-f(i);(beta*cR(1)-b1)*f(i)+(b2-beta*cR(2))*mR*g];
    for k = 1:length(ay)
        for j =1:3
            y(k,i,j) = (-D(j,k)*T(3,i)-A(j,k)*(T(1,i)+m*g)+C(j,k)*T(2,i))/...
                       (B(j,k)*(T(1,i)+m*g));
        end
        LT(k,i) = y(k,i,1) -max(y(k,i,2),y(k,i,3));
    end
end

contourf(f,ay,LT);

xlabel('f001');
ylabel('ay002');
zlabel('LT003');
hold on;
plot(f,(cE(1)*T(2,:)./(T(1,:)+m*g))+cE(2)+(T(3,:)./(T(1,:)+m*g)));
% 
% atest = linspace(-1,1,50)
% ctest = (atest.^2+2*L.*atest)./(2*L+atest);
% figure;
% plot(atest,ctest);
% 
% 
% 
% 
% 
% 
syms ay cy L l cx x y  mg fx fy Me real

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


k = simplify((D3*((ay-cy)*fx-cx*fy)+A3*fx+C3*fy)/(B3*fx));

sol2 = simplify(L+ay-((D2*Me+A2*fx+C2*fy)/(B2*fx)));
sol3 = simplify(L+ay-((D3*Me+A3*fx+C3*fy)/(B3*fx)));