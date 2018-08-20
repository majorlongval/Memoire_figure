function [minabs,minabs_i]= min_abs(v)
%MIN_ABS returns the value closest to zero in the vector v
minabs = inf;
minabs_i = 1;
for i =1:length(v)
    if abs(v(i)) < abs(minabs)
        minabs = v(i);
        minabs_i = i;
    end
end

