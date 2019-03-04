function [S_WFW] = Calc_S_WFW(Tx,Ty,Me,m,L,linput,ay,cy,cx)
%Calc_S_WFW calcul la surface du WFW
%   Detailed explanation goes here
ws = gen_wrench_set(Tx,Ty,Me);
for i1 =1:length(ws(:,1))
        LFM(:,:,i1) = [calc_lines2DDL(L,linput,ay,cy,cx,ws(i1,1),ws(i1,2),ws(i1,3),m);...
                      -1 0 0];

end


% 3.5) Gen a list of all lines
LIST_LINES = [];
for i2 =1:length(ws(:,1))
    LIST_LINES = [LIST_LINES;LFM(:,:,i2)];
end

% 4) Calculating all the intersection points 
intersection = cell(length(ws(:,1)),1);
for i3 =1:length(ws(:,1))-1
    for j3 =i3+1:length(ws(:,1))
   g_intersection =  all_inter_g(LFM(:,:,i3),LFM(:,:,j3));
   temp_i = intersection{i3,1};
   temp_i = [temp_i;g_intersection];
   intersection{i3,1} = temp_i;
   temp_j = intersection{j3,1};
   temp_j = [temp_j;g_intersection];
   intersection{j3,1} = temp_j;
    end
end

% 5) Filtering out bad data
 for i4 =1:length(intersection)
     temp = intersection{4};
     if ~isempty(temp)
     temp = temp(~isnan(temp(:,1)),:);
     temp = temp(~isinf(temp(:,1)),:);
     intersection{i4} = temp;
     end
 end

% 6) Ordering in lexicographic order the intersections and finding the
% vertices
vertices = [];
roundn = @(x,n) round(x.*10.^n)./10.^n;
for i5 =1:length(intersection)
    temp = intersection{i5};
    if ~isempty(temp)
    [Xp, Yp] = vec_acw_order_lin(temp(:,1),temp(:,2));
    for j5 =1:length(Xp)
        for k5 =1:length(LIST_LINES(:,1))
            bool(k5) = roundn(LIST_LINES(k5,1)*Xp(j5)+...
                             LIST_LINES(k5,2)*Yp(j5)+LIST_LINES(k5,3),8)<=0 ...
                             && Xp(j5) >= 0;
        end
        if all(bool== 1)
            vertices = [vertices;Xp(j5) Yp(j5)];
        end
    end
    end
end
if ~isempty(vertices)
        vertices = roundn(vertices,6);
    vertices = unique(vertices,'rows');
    vertices = sortrows(vertices,1);
    if length(vertices(:,1))<=2
    S_WFW = inf;
    else

[K,S_WFW] = convhull(vertices(:,1),vertices(:,2));
    end
else
    vertices = [];
    S_WFW = 0;
end






%%%%%%%%%%%%   SELF FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
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
for k = 1:length(Line_fact1(:,1))
    for l =1:length(Line_fact2(:,2))
        g_intersection = [g_intersection;calc_intersection(Line_fact1(k,:),Line_fact2(l,:))];
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
end

