function return_bool = ...
    check_COND_ETST(P1,P2,alpha, R, r, r_c, phi_c, h_c,tx, ty, tz, m, Mx,...
    My, Mz,z0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g =9.81;
Tz = tz + m*g;

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, r_c, phi_c, h_c);


roundn = @(x,n) round(x.*10.^n)./10.^n;
for i = 1:6
    A(i) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
    B(i) = tx*ay(i) + Tz*cy(i) + My*ey(i);
    C(i) = tx*(az(i)*z0+ak(i))+ty*(bz(i)*z0+bk(i))+Tz*ck(i)+...
        Mx*dk(i)+My*ek(i)+Mz*(hz(i)*z0+hk(i));
    if (roundn(A(i)*P1(1)+B(i)*P1(2)+C(i),4)>=0) &&...
            (roundn(A(i)*P2(1)+B(i)*P2(2)+C(i),4)>=0)
        list_elem(i) = 1;
    else
        list_elem(i) = 0;
    end
end
if all(list_elem(:) == 1)
    return_bool = 1;
else
    return_bool = 0;
end
end

