function minabs = min_abs(a,b)
%min_abs returns the closest value to zero.
%   a and b are the arguments 
if abs(a)<=abs(b)
    minabs = a;
else
    minabs = b;
end
end

