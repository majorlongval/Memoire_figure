function [ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
                calc_coeff(alpha, R, r, r_c, phi_c, h_c)
coeffX = ...
[ 0,-h_c,r_c*sin(phi_c)+2*r*cos(alpha)^2,1,0,0;...
  0,h_c, 2*r*cos(alpha)^2-r_c*sin(phi_c),-1,0,0;...
  0,-h_c,r_c*sin(phi_c)-r*cos(alpha)^2-(3^(1/2)*r*sin(2*alpha))/2,1,0,0;...
  0,h_c,-r_c*sin(phi_c)-r*cos(alpha)^2-(3^(1/2)*r*sin(2*alpha))/2,-1,0,0;...
  0,-h_c,r_c*sin(phi_c)-r*cos(alpha)^2+(3^(1/2)*r*sin(2*alpha))/2,1,0,0;...
  0,h_c,(3^(1/2)*r*sin(2*alpha))/2-r*cos(alpha)^2-r_c*sin(phi_c),-1,0,0];
coeffY = ...
[ h_c,0,r*sin(2*alpha)-r_c*cos(phi_c),0,1,0;...
 -h_c,0,r*sin(2*alpha)+r_c*cos(phi_c),0,-1,0;...
  h_c,0,3^(1/2)*r*cos(alpha)^2-r_c*cos(phi_c)-(r*sin(2*alpha))/2,0,1,0;...
 -h_c,0,r_c*cos(phi_c)-(r*sin(2*alpha))/2+3^(1/2)*r*cos(alpha)^2,0,-1,0;...
  h_c,0,-(r*sin(2*alpha))/2-r_c*cos(phi_c)-3^(1/2)*r*cos(alpha)^2,0,1,0;...
 -h_c,0,r_c*cos(phi_c)-(r*sin(2*alpha))/2-3^(1/2)*r*cos(alpha)^2,0,-1,0];
      
coeffZ = ...
[2*r*sin(alpha)^2-2*r-r_c*sin(phi_c),r_c*cos(phi_c)-r*sin(2*alpha),0,0,0,1;...
 2*r*sin(alpha)^2-2*r+r_c*sin(phi_c),-r*sin(2*alpha)-r_c*cos(phi_c),0,0,0,-1;...
 r*cos(alpha)^2-r_c*sin(phi_c)+(3^(1/2)*r*sin(2*alpha))/2,(r*sin(2*alpha))/2+r_c*cos(phi_c)-3^(1/2)*r*cos(alpha)^2,0,0,0,1;...
 r_c*sin(phi_c)+r*cos(alpha)^2+(3^(1/2)*r*sin(2*alpha))/2,(r*sin(2*alpha))/2-r_c*cos(phi_c)-3^(1/2)*r*cos(alpha)^2,0,0,0,-1;...
 r*cos(alpha)^2-r_c*sin(phi_c)-(3^(1/2)*r*sin(2*alpha))/2,(r*sin(2*alpha))/2+r_c*cos(phi_c)+3^(1/2)*r*cos(alpha)^2,0,0,0, 1;...
 r_c*sin(phi_c)+r*cos(alpha)^2-(3^(1/2)*r*sin(2*alpha))/2,(r*sin(2*alpha))/2-r_c*cos(phi_c)+3^(1/2)*r*cos(alpha)^2,0,0,0,-1];

coeffK = ...
[0, -2*R*h_c*cos(alpha),   R*r*cos(alpha) + 2*R*r_c*cos(alpha)*sin(phi_c),  2*R*cos(alpha),0, 0;...
0,  2*R*h_c*cos(alpha),   R*r*cos(alpha) - 2*R*r_c*cos(alpha)*sin(phi_c), -2*R*cos(alpha),0, 0;...
3^(1/2)*R*h_c*cos(alpha), R*h_c*cos(alpha), R*r*cos(alpha) - R*r_c*cos(alpha)*sin(phi_c) - 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c),-R*cos(alpha),3^(1/2)*R*cos(alpha), 0;...
-3^(1/2)*R*h_c*cos(alpha),-R*h_c*cos(alpha), R*r*cos(alpha) + R*r_c*cos(alpha)*sin(phi_c) + 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c),R*cos(alpha),-3^(1/2)*R*cos(alpha), 0;...
-3^(1/2)*R*h_c*cos(alpha),R*h_c*cos(alpha), R*r*cos(alpha) - R*r_c*cos(alpha)*sin(phi_c) + 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c),-R*cos(alpha),-3^(1/2)*R*cos(alpha), 0;...
3^(1/2)*R*h_c*cos(alpha), -R*h_c*cos(alpha), R*r*cos(alpha) + R*r_c*cos(alpha)*sin(phi_c) - 3^(1/2)*R*r_c*cos(alpha)*cos(phi_c),R*cos(alpha),3^(1/2)*R*cos(alpha), 0];

ay(:,1) = coeffY(:,1);
az(:,1) = coeffZ(:,1);
ak(:,1) = coeffK(:,1);
bx(:,1) = coeffX(:,2);
bz(:,1) = coeffZ(:,2);
bk(:,1) = coeffK(:,2);
cx(:,1) = coeffX(:,3);
cy(:,1) = coeffY(:,3);
ck(:,1) = coeffK(:,3);
dx(:,1) = coeffX(:,4);
dk(:,1) = coeffK(:,4);
ey(:,1) = coeffY(:,5);
ek(:,1) = coeffK(:,5);
hz(:,1) = coeffZ(:,6);
hk(:,1) = coeffK(:,6);






end
