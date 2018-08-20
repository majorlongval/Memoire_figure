function integ = ystat_surf_integ(xvect,yvect)
%YSTAT_surf_integ takes the integral of the surface of the ws from x0 to
%xend.
%   Detailed explanation goes here
integ = trapz(xvect,yvect-yvect(end));
end

