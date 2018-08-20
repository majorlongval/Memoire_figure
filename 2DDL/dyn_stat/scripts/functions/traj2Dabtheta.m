%% Traj simp

function [pos, acc, t,phix,phiy,rx,ry] = traj2Dabtheta(pos_ini, a, b, theta,phi, omega,T)
% TRAJ2DRXRYTHETA returns position and acceleation vectors of a trajectory
% defined by pos_ini, rx,ry, thera and omega. The transition time is the minimum
% time for which no tension is negative. The trajectory is planed as
% follows: transition to ellipse, three full turns and transition back to
% rest.


%% Calculating the values for phix and phiy

phix = atan2(b*sin(theta)*sin(phi)-a*cos(theta)*cos(phi),...
    a*cos(theta)*sin(phi)+b*sin(theta)*cos(phi));
phiy = atan2(a*sin(theta)*cos(phi)+b*cos(theta)*sin(phi),...
    b*cos(theta)*cos(phi)-a*sin(theta)*sin(phi));


rx = sqrt((a^2*cos(theta)^2+b^2*sin(theta)^2));
ry = sqrt((b^2*cos(theta)^2+a^2*sin(theta)^2));
%% Discretising time
tau = linspace(0,1,1000);
t1 = tau*T;
t2 = linspace(0,(6*pi),1000);
t3 = t1;
t  = [t1,t1(end)+t2,(t1(end)+t2(end)+t3)]; 
% setting the initial transition
delta1 = (6*tau.^5-15*tau.^4+10*tau.^3);
delta2 = (30*tau.^4-60*tau.^3+30*tau.^2)/T;
delta3 = (120*tau.^3-180*tau.^2+60*tau)/(T^2);
delta1y = ry*delta1;
delta1x = rx*delta1;
delta2x = rx*delta2;
delta2y = ry*delta2;
delta3x = delta3*rx;
delta3y = delta3*ry;
% setting the full trajectory
delta1y = [delta1y,ry*ones(1,length(t2))];
delta1x = [delta1x,rx*ones(1,length(t2))];
delta2x = [delta2x,zeros(1,length(t2))];
delta2y = [delta2y,zeros(1,length(t2))];
delta3x = [delta3x,zeros(1,length(t2))];
delta3y = [delta3y,zeros(1,length(t2))];
% setting the inverse transition trajectory
delta1y = [delta1y,ry*fliplr(delta1)];
delta1x = [delta1x,rx*fliplr(delta1)];
delta2x = [delta2x,rx*fliplr(delta2)];
delta2y = [delta2y,ry*fliplr(delta2)];
delta3x = [delta3x,rx*fliplr(delta3)];
delta3y = [delta3y,ry*fliplr(delta3)];


pos(1:2,:)   = pos_ini+ [delta1x.*sin(omega*t+phix);...
    delta1y.*sin(omega*t+phiy)];
acc(1:2,:) = [(delta3x-delta1x.*omega^2).*sin(omega*t+phix)+...
    2*delta2x.*omega.*cos(omega*t+phix);...
    (delta3y-delta1y.*omega^2).*sin(omega*t+phiy)+...
    2*delta2y.*omega.*cos(omega*t+phiy)];

% if we know the value of rx and ry and theta, a and b can be found with
% the following equations

% a^2 = (rx^2*cos(theta)^2-ry^s*sin(theta)^2)/(1-2*sin(theta));
% b^2 = (ry^2*cos(theta)^2-rx^s*sin(theta)^2)/(1-2*sin(theta));


end