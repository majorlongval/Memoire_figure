function p_inter = Calc_intersection(D1, D2)
%CALC_INTERSECTION calculates all the intersections between D1 and D2
%   D1 and D2 and lines that are given by Di = [Ai, Bi, Ci] such that
%   Aix+Biy+Ci = 0;
p_inter = (1/(D1(2)*D2(1)-D1(1)*D2(2)))*...
          [D1(3)*D2(2)-D1(2)*D2(3), D1(1)*D2(3)-D2(1)*D1(3)];
end

