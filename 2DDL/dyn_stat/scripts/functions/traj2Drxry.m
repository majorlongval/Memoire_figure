%% Traj simp

function [pos, acc, t] = traj2Drxry(pos_ini, rx, ry, omega, phix, phiy,tend)
% TRAJ2DRXRYTHETA returns position and acceleation vectors of a trajectory
% defined by pos_ini, rx,ry, thera and omega. The transition time is the minimum
% time for which no tension is negative. The trajectory is planed as
% follows: transition to ellipse, three full turns and transition back to
% rest.
%% Discretising time
t = linspace(tend,tend+7*pi/omega,1000);


% Calculating the position vector
for i =1:length(t)
    pos(:,i)   = pos_ini+ [rx*sin(omega*t(i)+phix);...
        ry*sin(omega*t(i)+phiy)];
    acc(:,i) = -omega^2*[rx*sin(omega*t(i)+phix);...
        ry*sin(omega*t(i)+phiy)];
end



end