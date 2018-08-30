function [ymin,ymax] = y_min_stat_tors(geo,x,tors,m)
[ay, cy, cx, L, l, A, B, C] = unpack_geo(geo);
D(1) = 0;
D(2) = (2*L+ay);
D(3) = -D(2);
g = 9.81;
q = tors(2)/(tors(1)+m*g);
for i =1:3
r(i) = (C(i)*tors(2)-D(i)*tors(3)-A(i)*(tors(1)+m*g))/(B(i)*(tors(1)+m*g));
y(i,:) = q*x(:)+r(i);
end

if (y(2,1)>y(3,1))
    ymin = y(2,:);
else
    ymin = y(3,:);
end

ymax = y(1,:);

    