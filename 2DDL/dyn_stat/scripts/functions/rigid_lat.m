function tau_y = rigid_lat(pos,geo)
% RIGID_LAT takes the position and the geometry and returns the latteral
% stifness
%   pos is the position of the end-effector, geo is a struct containing the
%   geometric parameters
[ay, cy, cx, L, l, A, B, C] = unpack_geo(geo);

for i =1:3
    tau_yt(i) = -9.81*(A(i)+B(i)*pos(2))/(C(i)+B(i)*pos(1));
end
tau_y_s = sort(tau_yt);
tau_y = tau_y_s(3)-tau_y_s(2);
