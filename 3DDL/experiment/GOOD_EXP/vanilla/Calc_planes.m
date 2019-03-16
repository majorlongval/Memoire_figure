function Planes = Calc_planes(alpha, R, r, r_c, phi_c, ...
                                         h_c,m,wrench_comb)
%Calc_planes calculates the equations of planes as a function of the geometric parameters
%in the form Ax+By+Cz+D = 0
g =9.81;
tx = wrench_comb(1);
ty = wrench_comb(2);
tz = wrench_comb(3);
Mx = wrench_comb(4);
My = wrench_comb(5);
Mz = wrench_comb(6);

Tz = tz + m*g;

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, r_c, phi_c, h_c);
for i = 1:6
    A(i,1) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
    B(i,1) = tx*ay(i) + Tz*cy(i) + My*ey(i);
    C(i,1) = tx*az(i)+ty*bz(i)+Mz*hz(i);
    D(i,1) = tx*ak(i)+ty*bk(i)+Tz*ck(i)+Mx*dk(i)+My*ek(i)+Mz*hk(i);

end
Planes = [A, B, C, D];
end

