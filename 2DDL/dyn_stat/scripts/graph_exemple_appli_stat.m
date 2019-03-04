%% Script to get the graph of tension in time with poly 5-5
clear all; close all; clc;
addpath('functions')
g = 9.81;

L = 5;
l = 0.8;

cx = 0.2;
cy =-0.1;
ay = 0.3;

mR = 2;
mC = 1;

dx = 0.4;
dy = -0.2;
ex = 0.3;
ey = -0.3;

geo = pack_geo(ay,cy,cx,5,1);
[ay, cy, cx, L, l,A, B, C] = unpack_geo(geo);
D(1) = 0;
D(2) = 2*L+ay;
D(3) = -D(2);

% Points 
PI = [3;0.5];
PCp = [4;0];
PC = [5;0];
PE = [2.5;1];

TorsA = [0;0;0];
% Calculating the time from PI to PCp
T1 = T_over(TorsA,A,B,C,D,mR,PI,PCp);
% Calculating the time from PCp to PC
T2 = T_over(TorsA,A,B,C,D,mR,PCp,PC);
TorsB = mC*g*[1;0;dy];
% Calculating the time from PC to PCp
T3 = T_over(TorsB,A,B,C,D,mR,PC,PCp);
% Calculating the time from PCp to PE
T4 = T_over(TorsB,A,B,C,D,mR,PCp,PE);


% Calculating the trajectories and the tensions
%1)
[pos1, acc1, t1] = calc_pos_acc(PI,PCp,T1);
for i =1:length(pos1(1,:))
    tens1(:,i) = tension(mR,geo,pos1(:,i),acc1(:,i),TorsA);
end
%2)
[pos2, acc2, t2] = calc_pos_acc(PCp,PC,T2);
for i =1:length(pos1(1,:))
    tens2(:,i) = tension(mR,geo,pos2(:,i),acc2(:,i),TorsA);
end
%3)
[pos3, acc3, t3] = calc_pos_acc(PC,PCp,T3);
for i =1:length(pos1(1,:))
    tens3(:,i) = tension(mR,geo,pos3(:,i),acc3(:,i),TorsB);
end
%4)
[pos4, acc4, t4] = calc_pos_acc(PCp,PE,T4);
for i =1:length(pos1(1,:))
    tens4(:,i) = tension(mR,geo,pos4(:,i),acc4(:,i),TorsB);
end



% combining the things
t2 = t2 + t1(end);
t3 = t3 + t2(end);
t4 = t4 + t3(end);
t = [t1,t2,t3,t4];

tens = [tens1,tens2,tens3,tens4];
p1 = plot(t,tens(1,:),'-b');
hold on;
p2 = plot(t,tens(2,:),'-r');
p3 = plot(t,tens(3,:),'-k');
grid on;
ph1 = plot([T1 T1],[0 35],':k');
ph1.LineWidth = 2;
ph2 = plot([T1+T2 T1+T2],[0 35],':k');
ph2.LineWidth = 2;
ph3 = plot([T1+T2+T3 T1+T2+T3],[0 35],':k');
ph3.LineWidth = 2;

saveas(gca,'time_tens_exemple.svg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [q1,q3] = calc_poly_5(T)
t= linspace(0,T,1000);
tau = t./T;
q1 = 6*tau.^5 -15*tau.^4 + 10*tau.^3;
q3 = (120*tau.^3 -180*tau.^2 + 60*tau)/(T^2);
end

function [pos, acc, t] = calc_pos_acc(pI,pF,T)
[q1,q3] = calc_poly_5(T);
t = linspace(0,T,length(q1));
pos = pI+q1.*(pF-pI);
acc = q3.*(pF-pI);
end


function allbool = test_cond_T(Tors,A,B,C,D,m,pI,pF,T)
[q1,q3] = calc_poly_5(T);
xI = pI(1);
xF = pF(1);
yI = pI(2);
yF = pF(2);
dy = yF-yI;
dx = xF-xI;
Tx = Tors(1) + m*9.81;
ty = Tors(2);
Me = Tors(3);
for i =1:length(q1)
    for j =1:3
        u(i,j) = (A(j)+B(j)*(yI+q1(i)*dy))*Tx-...
                 (C(j)+B(j)*(xI+q1(i)*dx))*ty+D(j)*Me;
        v(i,j) = m*(C(j)*dy-A(j)*dx+B(j)*(xI*yF-xF*yI))*q3(i);
    end
   if u(i,1)+v(i,1)<0 && u(i,2)+v(i,2)<0 && u(i,3)+v(i,3)<0
       bool(i) = 1;
   else
       bool(i) = 0;
   end  
end
allbool = all(bool);
end

function T_opt = min_T(Tors,A,B,C,D,m,pI,pF)
T = 0.001;
allbool = 0;
while allbool== 0
    T = T+0.001;
    allbool = test_cond_T(Tors,A,B,C,D,m,pI,pF,T);
end
T_opt = T;
end

function Tover = T_over(Tors,A,B,C,D,m,pI,pF)
xI = pI(1);
xF = pF(1);
yI = pI(2);
yF = pF(2);
dy = yF-yI;
dx = xF-xI;
Tx = Tors(1) + m*9.81;
ty = Tors(2);
Me = Tors(3);
for i =1:3
    v(i) = abs(m*(C(i)*dy-A(i)*dx+B(i)*(xI*yF-xF*yI)))*10*sqrt(3)/3;
    uI(i) = (A(i)+B(i)*yI)*Tx-(C(i)+B(i)*xI)+D(i)*Me;
    uF(i) = (A(i)+B(i)*yF)*Tx-(C(i)+B(i)*xF)+D(i)*Me;
    T1(i) = sqrt(-v(i)/uI(i));
    T2(i) = sqrt(-v(i)/uF(i));
end
Tover = max([T1(1),T1(2),T1(3),T2(1),T2(2),T2(3)]);



end