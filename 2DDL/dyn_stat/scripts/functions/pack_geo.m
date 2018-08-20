function geo = pack_geo(ay,cy,cx,L,l)
%PACK_GEO packs the geometric parameters in a struct
%   Detailed explanation goes here
geo.ay = ay;
geo.cy = cy;
geo.cx = cx;
geo.L = L;
geo.l = l;


geo.A(1) = -L;
geo.A(2) = ay*(cy-L)+2*L*cy-l*(ay+L);
geo.A(3) = -ay*(cy-L)-2*L*cy-l*(ay+L);

geo.B(1) = 1;
geo.B(2) = ay-l;
geo.B(3) = -(ay+l);

geo.C(1) = 0;
geo.C(2) = cx*(2*L+ay);
geo.C(3) = -geo.C(2);
end

