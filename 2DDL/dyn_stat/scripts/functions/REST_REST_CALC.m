%% Calculating the transition period for point to point inside ETST
clear all; close all; clc;

%addpath('functions');


m = 1; 
g = 9.81;

L = 5;
l = 0.8;

cx = 0.2;
cy =0.1;
ay = -0.1;

tx = 0;
ty = 0;
tz = 0;


[ay, cy, cx, L, l,A, B, C] = unpack_geo(pack_geo(ay,cy,cx,5,1));



Tors = [0 0 0];




D(1) = 0;
D(2) = 2*L+ay;
D(3) = -D(2);

% pente = Tors(2)/(m*g+Tors(1));
% 
% b(1:3) = (C(1:3).*Tors(2)-A(1:3).*(m*g+Tors(1))-D(1:3).*Tors(3))./...
%          (B(1:3).*(m*g+Tors(1)));
% 
% x = linspace(0, 5,1000);
% 
% y(:,1) = pente*x+b(1);
% y(:,2) = pente*x+b(2);
% y(:,3) = pente*x+b(3);
% 
% y2(:,1) = ones(1,1000)*-A(1)/B(1);
% y2(:,2) = ones(1,1000)*-A(2)/B(2);
% y2(:,3) = ones(1,1000)*-A(3)/B(3);
% 
% figure;
% 
% ph = patch([x(1) x(end) x(end) x(1)],[y(1,2) y(end,2) y(end,1) y(1,1)],...
%       [0.8,0.8,0.8],'EdgeColor','none');
% 
% ph.FaceAlpha = 0.6;
% % p0 = fill([x,x(end),fliplr(x),x(1)],[y(:,3)',...
% %     y(end,1),fliplr(y(:,1)'),y(1,1)],[0.8,0.8,0.8],'EdgeColor','none');
% hold on;
% p1 = plot(x,y(:,1),'-b','LineWidth',2);
% p2 = plot(x,y(:,2),'-r','LineWidth',2);
% p3 = plot(x,y(:,3),'-k','LineWidth',2);
% 
% 
% %patch([x(1) x(end) x(end) x(1)],[y(1,2) y(end,2) y(end,1) y(1,1)],'r');
% 
% camroll(-90);
% grid on;
% axis([0 5 -10 10]);
% ax=gca;
% ax.YAxisLocation = 'Right';
% xlabel('yyy');
% ylabel('xxx');
% title('TITLE');
% 
% 
% Plotting two points 
pI = [1;-2];
pF = [4.5;4];
% 
% plot(pI(1),pI(2),'*r');
% hold on;
% plot(pF(1),pF(2),'*b');
% plot([pI(1) pF(1)],[pI(2) pF(2)],'--k');
% 
% 
% saveas(gca,'traj_rest_rest.svg');
% Calculating the trajectory
% 
% %method 1: Numerical itteration
% T_opt = min_T(Tors,A,B,C,D,m,pI,pF);

%method 2: overestimation method
Tover = T_over(Tors,A,B,C,D,m,pI,pF);


% Calculating the tensions for many cases 
% geo = pack_geo(ay,cy,cx,L,l);
% [pos_opt, acc_opt, t_opt] = calc_pos_acc(pI,pF,T_opt);
% for i =1:length(pos_opt(1,:))
%     tens_opt(:,i) = tension(m,geo,pos_opt(:,i),acc_opt(:,i),Tors');
% end
% figure;
% plot(t_opt,tens_opt(1,:),'-r');
% hold on;
% plot(t_opt,tens_opt(2,:),'-b');
% plot(t_opt,tens_opt(3,:),'-k');
% grid on;

[pos_over, acc_over, t_over] = calc_pos_acc(pI,pF,Tover);
for i =1:length(pos_opt(1,:))
    tens_over(:,i) = tension(m,geo,pos_over(:,i),acc_over(:,i),Tors');
end
plot(t_over,tens_over(1,:),'--r');
plot(t_over,tens_over(2,:),'--b');
plot(t_over,tens_over(3,:),'--k');

axis([0 4 0 25]);

plot([Tover Tover],[0 25],':k');
plot([T_opt T_opt],[0 25],':k');
xlabel('xxx'); ylabel('yyy');

legend('legend1','legend2','legend3','legend4','legend5','legend6')

saveas(gca,'traj_rest_rest_tens.svg');

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