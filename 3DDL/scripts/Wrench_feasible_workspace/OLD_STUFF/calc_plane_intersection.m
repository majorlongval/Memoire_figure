function X = calc_plane_intersection(p1,p2,p3)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
M =[p1(1),p1(2),p1(3);...
    p2(1),p2(2),p2(3);...
    p3(1),p3(2),p3(3)];
if rank(M) == 3
    d = -[p1(4);p2(4);p3(4)];
    X = M\d;
else
    X = [];
    
    
    
end

