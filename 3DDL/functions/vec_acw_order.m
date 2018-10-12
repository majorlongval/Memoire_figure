function [Xp, Yp] = vec_acw_order(X,Y)
%VEC_CW_ORDER puts the vectors in an anti clockwise order so starting at
%the vector with the smallest angle from the x axis.
%   X is a vector containing the x coordinates of the vectors
%   Y is a vector containing the y coordinates of the vectors
%   Xp is a vector containing the x coordinates of the vector in a acw
%   order.
%   Yp is a vector containing the y coordinates of the vector in a acw
%   order.

X_bary  = mean(X);
Y_bary  = mean(Y);

X_comp = X-X_bary;
Y_comp = Y-Y_bary;

for i =1:length(X)
    mat(i,:) = [X(i),Y(i),atan2(Y_comp(i),X_comp(i))];
end

B = sortrows(mat,3,'descend');

Xp = B(:,1); Yp = B(:,2);
end

