%% Script to plot the WFW in the cartesian space.
clear all; close all; clc;
m =1;
% 1) Generating the wrench set 
Tx = [-0.06,0.06]*9.81*m;%[-0.1 0.1];
Ty = [-0.05,0.05]*9.81*m;
Me = [-0.05, 0.05]*9.81*m;
ws = gen_wrench_set(Tx,Ty,Me);

% 2) Setting the geometric values 
m = 1;
L = 2; 
l = 1;
ay = 3;
cy = ay/2;
cx = 0.12;

% 3) Gen the Line_fact 
for i =1:length(ws(:,1))
        LFM(:,:,i) = [calc_lines2DDL(L,l,ay,cy,cx,ws(i,1),ws(i,2),ws(i,3),m);...
                      -1 0 0];

end


% 3.5) Gen a list of all lines
LIST_LINES = [];
for i =1:length(ws(:,1))
    LIST_LINES = [LIST_LINES;LFM(:,:,i)];
end

% 4) Calculating all the intersection points 
intersection = cell(length(ws(:,1)),1);
for i =1:length(ws(:,1))-1
    for j =i+1:length(ws(:,1))
   g_intersection =  all_inter_g(LFM(:,:,i),LFM(:,:,j));
   temp_i = intersection{i};
   temp_i = [temp_i;g_intersection];
   intersection{i} = temp_i;
   temp_j = intersection{j};
   temp_j = [temp_j;g_intersection];
   intersection{j} = temp_j;
    end
end

% 5) Filtering out bad data
 for i =1:length(intersection)
     temp = intersection{i};
     temp = temp(~isnan(temp(:,1)),:);
     temp = temp(~isinf(temp(:,1)),:);
     intersection{i} = temp;
 end

% 6) Ordering in lexicographic order the intersections and finding the
% vertices
vertices = [];
roundn = @(x,n) round(x.*10.^n)./10.^n;
for i =1:length(intersection)
    temp = intersection{i};
    [Xp, Yp] = vec_acw_order_lin(temp(:,1),temp(:,2));
    for j =1:length(Xp)
        for k =1:length(LIST_LINES(:,1))
            bool(k) = roundn(LIST_LINES(k,1)*Xp(j)+...
                             LIST_LINES(k,2)*Yp(j)+LIST_LINES(k,3),5)<=0 ...
                             && Xp(j) >= 0;
        end
        if all(bool== 1)
            vertices = [vertices;Xp(j) Yp(j)];
        end
    end
end
vertices = roundn(vertices,3);
vertices = unique(vertices,'rows');
 vertices = sortrows(vertices,1);


% 7) Plotting lines
x = 0:0.01:5;
for i =1:length(LIST_LINES(:,1))
    for j =1:length(x)
        y(j,i) = -(LIST_LINES(i,1)*x(j)+LIST_LINES(i,3))/LIST_LINES(i,2);
    end
    plot(x,y(:,i),'-b');
    hold on;
end

plot(vertices(:,1),vertices(:,2),'*g');
% ph = patch(vertices(:,1),vertices(:,2),'r');
% ph.FaceAlpha = 0.5;
grid on;
camroll(-90);
xlabel('XXX');
ylabel('YYY');

axis([0 5 -6 4]);
% for i =1:length(intersection)
%     temp=intersection{i}
%     plot(temp(:,1),temp(:,2),'*r');
% end

saveas(gca,'all_lines_2DDL.svg');
[k,S] = convhull(vertices(:,1),vertices(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Line_fact = calc_lines2DDL(L,l,ay,cy,cx,tx,ty,Me,m)
g = 9.81; % m/s^²
A1 = -L;
A2 = ay*(cy-L)+2*L*cy-l*(ay+L);
A3 = ay*(L-cy)-2*L*cy-l*(ay+L);

B1 = 1;
B2 = ay-l;
B3 = -(ay+l);

C1 = 0;
C2 = cx*(2*L+ay);
C3 = -C2;

D1 = 0;
D2 = 2*L+ay;
D3 = -D2;



a1 = -B1*ty;
a2 = -B2*ty;
a3 = -B3*ty;
b1 = B1*(m*g+tx);
b2 = B2*(m*g+tx);
b3 = B3*(m*g+tx);
c1 = A1*(m*g+tx)-C1*ty+D1*Me;
c2 = A2*(m*g+tx)-C2*ty+D2*Me;
c3 = A3*(m*g+tx)-C3*ty+D3*Me;


% q = ty/(tx+m*g);
% r1 = (C1*ty-D1*Me-A1*(tx+m*g))/(B1*(tx+m*g));
% r2 = (C2*ty-D2*Me-A2*(tx+m*g))/(B2*(tx+m*g));
% r3 = (C3*ty-D3*Me-A3*(tx+m*g))/(B3*(tx+m*g));

Line_fact = [a1, b1, c1;...
             a2, b2, c2;...
             a3, b3, c3];
end

function intersection = calc_intersection(D1,D2)
intersection = (1/(D1(2)*D2(1)-D1(1)*D2(2)))*...
          [D1(3)*D2(2)-D1(2)*D2(3), D1(1)*D2(3)-D2(1)*D1(3)];
end

function g_intersection = all_inter_g(Line_fact1,Line_fact2)
g_intersection = [];
for i = 1:length(Line_fact1(:,1))
    for j =1:length(Line_fact2(:,2))
        g_intersection = [g_intersection;calc_intersection(Line_fact1(i,:),Line_fact2(j,:))];
    end
end

end

function wrench_set = gen_wrench_set(Tx,Ty,Me)
wrench_set = [];
for i = 1:length(Tx)
    for j =1:length(Ty)
        for k=1:length(Me)
            wrench_set = [wrench_set;[Tx(i),Ty(j),Me(k)]];
        end
    end
end
end

function [Xp, Yp] = vec_acw_order_lin(X,Y)
%VEC_CW_ORDER puts the vectors in an clockwise order so starting at
%the vector with the smallest angle from the x axis.
%   X is a vector containing the x coordinates of the vectors
%   Y is a vector containing the y coordinates of the vectors
%   Xp is a vector containing the x coordinates of the vector in a acw
%   order.
%   Yp is a vector containing the y coordinates of the vector in a acw
%   order.

mat = [X,Y];
if all(X == X(1))
    mat_p = sortrows(mat,2);
else
    mat_p = sortrows(mat,1);
end
Xp = mat_p(:,1);
Yp = mat_p(:,2);
end