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

