function List_lines = Calc_lines(alpha, R, r, r_c, phi_c, ...
                                         h_c,m,z0,wrench_comb)%,rhomax,rhomin,fmin,fmax)
%CALC_LINES Calculates the A B and C coefficients of the lines associated
%with a set of wrench 
%   List_lines is a matrix of the form [A B C] where A B C are column
%   vectors of length 6.
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
    C(i,1) = tx*(az(i)*z0+ak(i))+ty*(bz(i)*z0+bk(i))+Tz*ck(i)+...
        Mx*dk(i)+My*ek(i)+Mz*(hz(i)*z0+hk(i));%-(6*R*r*z0*cos(alpha)*fmin/rhomin(i));
end
% for i = 7:12
%     A(i,1) = -A(i-6,:);
%     B(i,1) = -B(i-6,:);
%     C(i,1) = -(tx*(az(i-6)*z0+ak(i-6))+ty*(bz(i-6)*z0+bk(i-6))+Tz*ck(i-6)+...
%         Mx*dk(i-6)+My*ek(i-6)+Mz*(hz(i-6)*z0+hk(i-6)))+(6*R*r*z0*cos(alpha)*fmax/rhomax(i-6));
% end
List_lines = [A,B,C];
end

