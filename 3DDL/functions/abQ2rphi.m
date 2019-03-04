function [rx,ry,rz,phix,phiy,phiz] = abQ2rphi(a,b,Q,phi)
%abQ2rphi Takes a typical elliptic trajectory in 3d and gives the simpified
%form of it
%   a is the major axis
%   b is the minor axis
%   Q is the rotation matrix
%   phi is the offset
%   The input form must be of the form p =
%   Q[acos(theta+phi);bsin(theta+phi);0],where theta is the variable
q11 = Q(1,1); q12 = Q(1,2); q13 = Q(1,3); 
q21 = Q(2,1); q22 = Q(2,2); q23 = Q(2,3);
q31 = Q(3,1); q32 = Q(3,2); q33 = Q(3,3);
cp = cos(phi); sp = sin(phi);

cp = cos(phi); sp = sin(phi);

MM = [q11*a*cp+q12*b*sp q12*b*cp-q11*a*sp;...
      q21*a*cp+q22*b*sp q22*b*cp-q21*a*sp;...
      q31*a*cp+q32*b*sp q32*b*cp-q31*a*sp];

phix = atan2(MM(1,1),MM(1,2));
phiy = atan2(MM(2,1),MM(2,2));
phiz = atan2(MM(3,1),MM(3,2));

rx =  sqrt(MM(1,1)^2+MM(1,2)^2);
ry =  sqrt(MM(2,1)^2+MM(2,2)^2);
rz =  sqrt(MM(3,1)^2+MM(3,2)^2);

end

