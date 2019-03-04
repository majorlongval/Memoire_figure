%% Traj simp

function [pos, acc] = traj2Drxry(pos_ini, rx, ry, omega, phix, phiy)
% TRAJ2DRXRYTHETA returns position and acceleation vectors of a trajectory
% defined by pos_ini, rx,ry, thera and omega. The transition time is the minimum
% time for which no tension is negative. The trajectory is planed as
% follows: transition to ellipse, three full turns and transition back to
% rest.
           
rx  = 
%% Discretising time
t = linspace(0,2*pi,1000);


% Calculating the position vector
for i =1:length(t)
    pos(:,i)   = pos_ini+ [rx*sin(omega*t(i)+phix);...
        ry*sin(omega*t(i)+phiy)];
    acc(:,i) = -omega^2*[rx*sin(omega*t(i)+phix);...
        ry*sin(omega*t(i)+phiy)];
end

% if we know the value of rx and ry and theta, a and b can be found with
% the following equations

% a^2 = (rx^2*cos(theta)^2-ry^s*sin(theta)^2)/(1-2*sin(theta));
% b^2 = (ry^2*cos(theta)^2-rx^s*sin(theta)^2)/(1-2*sin(theta));


end