function return_bool = check_COND_INTER(P1,P2,alpha, R, r, r_c, ...
                                        phi_c, h_c,tx, ty, tz,...
                                        m, Mx,My, Mz,z0)
%CHECK_COND_INTER Checks if the point respects both the static worksapce
%and the static workspace with external wrench 
%   Detailed explanation goes here
return_bool_1 = ...
    check_COND_ETST(P1,P2,alpha, R, r, r_c, phi_c, h_c,tx, ty, tz, m, Mx,...
    My, Mz,z0);
return_bool_2 = check_cond_new(P1,P2,R,r,r_c,alpha,phi_c);

if return_bool_1 && return_bool_2
    return_bool = 1;
else
    return_bool = 0;
end

