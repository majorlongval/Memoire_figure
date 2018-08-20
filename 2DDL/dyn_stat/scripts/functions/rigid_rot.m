function tau_r = rigid_rot(pos,geo)
% RIGID_LAT takes the position and the geometry and returns the latteral
% stifness
%   pos is the position of the end-effector, geo is a struct containing the
%   geometric parameters
[ay, cy, cx, L, l, A, B, C] = unpack_geo(geo);

for i =1:2
    tau_rt(i) = abs(-9.81*(A(i+1)+B(i+1)*pos(2))/(2*L+ay));
end
tau_r = min(tau_rt);