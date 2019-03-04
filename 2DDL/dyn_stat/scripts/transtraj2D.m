function [posup,posdown, accup,accdown,tup,tdown] = transtraj2D(pos_ini, rx, ry,phix,phiy, omega, T)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tau = linspace(0,1,1000);
t = tau*T;
delta1 = (6*tau.^5-15*tau.^4+10*tau.^3);
delta1y = ry*delta1;
delta1x = rx*delta1;
delta2 = (30*tau.^4-60*tau.^3+30*tau.^2)/T;
delta2x = rx*delta2;
delta2y = ry*delta2;
delta3 = (120*tau.^3-180*tau.^2+60*tau)/(T^2);
delta3x = delta3*rx;
delta3y = delta3*ry;


    posup(1:2,:)   = pos_ini+ [delta1x.*sin(omega*t+phix);...
        delta1y.*sin(omega*t+phiy)];
    accup(1:2,:) = [(delta3x-delta1x.*omega^2).*sin(omega*t+phix)+...
        2*delta2x.*omega.*cos(omega*t+phix);...
        (delta3y-delta1y.*omega^2).*sin(omega*t+phiy)+...
        2*delta2y.*omega.*cos(omega*t+phiy)];
    
posdown = fliplr(posup);
accdown = fliplr(accup);
    
% delta1 = (6*tau.^5-15*tau.^4+10*tau.^3);
% delta1 = fliplr(delta1);
% delta1y = ry*delta1;
% delta1x = rx*delta1;
% delta2 = (30*tau.^4-60*tau.^3+30*tau.^2)/T;
% delta2 = fliplr(delta2);
% delta2x = rx*delta2;
% delta2y = ry*delta2;
% delta3 = (120*tau.^3-180*tau.^2+60*tau)/(T^2);
% delta3 = fliplr(delta3);
% delta3x = delta3*rx;
% delta3y = delta3*ry;
% 
% 
%     posdown(1:2,:)   = pos_ini+ [delta1x.*sin(omega*t+phix);...
%         delta1y.*sin(omega*t+phiy)];
%     accdown(1:2,:) = [(delta3x-delta1x.*omega^2).*sin(omega*t+phix)+...
%         2*delta2x.*omega.*cos(omega*t+phix);...
%         (delta3y-delta1y.*omega^2).*sin(omega*t+phiy)+...
%         2*delta2y.*omega.*cos(omega*t+phiy)];
%     
tup = t;
tdown = t;
end

