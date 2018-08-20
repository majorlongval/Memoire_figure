function tau = tension_stat(pos,gp)
%TENSION_STAT Calculates the static tension/unit mass in each cable as a function of
%the position and the geometric parameters of the robot
%   pos is the [x;y] vector of position of the end-effector
%   gp is a struct containing the geometric parameters of the robot

% Importing values
L = gp.L;
l = gp.l;
a1 = [gp.ax;gp.ay];
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
delta1    = (a1-c)'*E*e1;
delta2    = (a2-c)'*E*e2;
delta3    = (a3-c)'*E*e2;

M = [e1 e2 e2;delta1 delta2 delta3];
tau = M\gam;


end

