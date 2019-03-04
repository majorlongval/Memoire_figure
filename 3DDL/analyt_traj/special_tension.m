clear all; close all; clc;

syms ayi azi bxi bzi cki cyi cxi aki bki cki dxi dyi dzi zki bki cki dki eyi ...
     eki hki hzi ...
     xc yc zc tx ty tz Mx My Mz m g omega phix phiy phiz  kx kz  Tz real 
Lambda = [0 bxi cxi;ayi 0 cyi;azi bzi 0];
Ni = [Lambda,[dxi 0 0;0 eyi 0;0 0 hzi]];
lambda = [aki;bki;cki];
ni = [aki;bki;cki;dki;eki;hki];
t = [tx;ty;Tz;Mx;My;Mz];
pc = [xc;yc;zc];
k2 = [kx*sin(phix);sin(phiy);kz*sin(phiz)];
k1 = [kx*cos(phix);cos(phiy);kz*cos(phiz)];
ui = m*omega^2*(pc'*Lambda+lambda')+t'*Ni';
vi = (pc'*Ni+ni')*t;
Thetai = sqrt((ui*k1)^2 + (ui*k2)^2);

elem1 = subs(ui*k1,[ayi,aki,bxi,bki],[0,0,0,0]);
elem1 = subs(elem1,omega,sqrt((Tz)/(m*zc)));
elem1 = subs(elem1,[Mx,My,Mz],[0,0,0]);
elem1 = elem1^2;
elem2 = subs(ui*k2,[ayi,aki,bxi,bki],[0,0,0,0]);
elem2 = subs(elem2,omega,sqrt((Tz)/(m*zc)));
elem2 = subs(elem2,[Mx,My,Mz],[0,0,0]);
elem2 = elem2^2;
% elem9 = simplify(elem4+elem8);
elem3 = subs(vi,[ayi,aki,bxi,bki],[0,0,0,0]);
elem3 = subs(elem3,[Mx,My,Mz],[0,0,0]);
% elem12 = simplify(elem11/sqrt(elem9));