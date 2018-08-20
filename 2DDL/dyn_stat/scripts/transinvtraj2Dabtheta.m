function [pos, acc] = transtraj2Dabtheta(pos_ini, a, b, theta,phi, omega, T)


phix = atan2(b*sin(theta)*sin(phi)-a*cos(theta)*cos(phi),...
    a*cos(theta)*sin(phi)+b*sin(theta)*cos(phi));
phiy = atan2(a*sin(theta)*cos(phi)+b*cos(theta)*sin(phi),...
    b*cos(theta)*cos(phi)-a*sin(theta)*sin(phi));


rx = sqrt((a^2*cos(theta)^2+b^2*sin(theta)^2));
ry = sqrt((b^2*cos(theta)^2+a^2*sin(theta)^2));


tau = linspace(0,1,1000);
t = tau*T;
delta1 = (6*tau.^5-15*tau.^4+10*tau.^3);
delta1 = fliplr(delta1);
delta1y = ry*delta1;
delta1x = rx*delta1;
delta2 = (30*tau.^4-60*tau.^3+30*tau.^2)/T;
delta2 = fliplr(delta2);
delta2x = rx*delta2;
delta2y = ry*delta2;
delta3 = (120*tau.^3-180*tau.^2+60*tau)/(T^2);
delta3 = fliplr(delta3);
delta3x = delta3*rx;
delta3y = delta3*ry;


    pos(1:2,:)   = pos_ini+ [delta1x.*sin(omega*t+phix);...
        delta1y.*sin(omega*t+phiy)];
    acc(1:2,:) = [(delta3x-delta1x.*omega^2).*sin(omega*t+phix)+...
        2*delta2x.*omega.*cos(omega*t+phix);...
        (delta3y-delta1y.*omega^2).*sin(omega*t+phiy)+...
        2*delta2y.*omega.*cos(omega*t+phiy)];
end
