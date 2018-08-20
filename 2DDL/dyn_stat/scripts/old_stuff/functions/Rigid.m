function [tau_y,Msol] = Rigid(pos,gp)
%Rigid gives an indication of the laterral and rotational stifness of the mechanism
%   tau_y is the latteral stifness and M is the rotational stifness


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating the vectors and moments
L = gp.L;
l = gp.l;
a1 = [0;gp.ay];
a2 = [0;-gp.l];
a3 = [0;gp.l];
b1 = [0;-gp.L];
b2 = [0;gp.L-gp.l];
b3 = [0;gp.L+gp.l];
c  = [gp.cx;gp.cy];

g = 9.81;

E = [0,-1;1,0];
gam = [g;0;0];

if size(pos) == [1,2]
    pos = pos';
end

e1 = (pos+a1-b1)/norm(pos+a1-b1);
e2 = (pos+a2-b2)/norm(pos+a2-b2);
delta1    = (c-a1)'*E*e1;
delta2    = (c-a2)'*E*e2;
delta3    = (c-a3)'*E*e2;

% Latteral stifness
ey = [0;1];

M1 = [e2 e2 ey;delta2 delta3 0];
M2 = [e1 e2 ey;delta1 delta3 0];
M3 = [e1 e2 ey;delta1 delta2 0];

sol1 = M1\gam;
sol2 = M2\gam;
sol3 = M3\gam;

tau_y = min([abs(sol1(3)),abs(sol2(3)),abs(sol3(3))]);


% Rotational stifness

er = [0; 0];

% M1 = [e2 e2 er;delta2 delta3 1];
M2 = [e1 e2 er;delta1 delta3 1];
M3 = [e1 e2 er;delta1 delta2 1];

% sol1 = M1\gam;
sol2 = M2\gam;
sol3 = M3\gam;

Msol = min([abs(sol2(3)),abs(sol3(3))]);


end

