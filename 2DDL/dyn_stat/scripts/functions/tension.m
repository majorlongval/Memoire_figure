function tension = tension(m,geo,pos,acc,tors)
% TENSION returns the tensions in the cables given m, the mass, geo, the
% geometric parameters, pos, the position and acc the acceleration

g = 9.81; %m/s^2
% Unpacking geo
[ay, cy, cx, L, l,A, B, C] = unpack_geo(geo);

% Calculating the constant vectors
a1 = [0;ay];
a2 = [0;-l];
a3 = [0;l];

b1 = [0;-L];
b2 = [0;L-l];
b3 = [0;L+l];

c = [cx;cy];

E = [0,-1;1,0];


% calculating the Jacobian Matrix
e1 = (pos+a1-b1)/norm(pos+a1-b1);
e2 = (pos+a2-b2)/norm(pos+a2-b2);
delta1    = (a1-c)'*E*e1;
delta2    = (a2-c)'*E*e2;
delta3    = (a3-c)'*E*e2;
M         = [e1 e2 e2;delta1 delta2 delta3];

% Calculating the gamma vector
gam = tors+[m*(g-acc(1));-m*acc(2);0];

tension = M\gam;

end

