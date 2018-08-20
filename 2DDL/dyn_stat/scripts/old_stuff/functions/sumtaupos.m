function sumtau = sumtaupos(tau)
%sumtaupos return the sum of all the tensions if all the tensions are
%positive
if tau(1)>0 & tau(2)>0 & tau(3)>0
    sumtau = tau(1)+tau(2)+tau(3);
end

