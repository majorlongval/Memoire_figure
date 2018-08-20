function [ay, cy, cx, L, l,A, B, C] = unpack_geo(geo)
%UNPACK_GEO unpacks the structure to variables
ay = geo.ay;
cy = geo.cy;
cx = geo.cx;
L  = geo.L;
l  = geo.l;
A  = geo.A;
B  = geo.B;
C  = geo.C;
end

