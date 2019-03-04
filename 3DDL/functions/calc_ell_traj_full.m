function [t,pos,acc,it1,it2,it3] = calc_ell_traj_full(pc,omega,rx,ry,rz,px,py,pz,T,nbt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tau = linspace(0,1,1000);
t1 = tau*T;
t2 = t1(end)+linspace(0,(nbt*pi),1000); 
t3 = t2(end)+t1;
%t  = t2;%[t1,t1(end)+t2,(t1(end)+t2(end)+t3)];
% setting the initial transition
delta1 = (6*tau.^5-15*tau.^4+10*tau.^3);
delta2 = (30*tau.^4-60*tau.^3+30*tau.^2)/T;
delta3 = (120*tau.^3-180*tau.^2+60*tau)/(T^2);

Kx = rx/ry; Kz = rz/ry;
cpx = cos(px);
spx = sin(px);
cpy = cos(py);
spy = sin(py);
cpz = cos(pz);
spz = sin(pz);

K   = [Kx*cpx Kx*spx;cpy spy;Kz*cpz Kz*spz];
E  = [0 1;-1 0];

s1 = [sin(omega*t1);cos(omega*t1)];
pos1 = pc+ry*[delta1;delta1;delta1].*(K*s1);
acc1 = ry*([(delta3-omega^2*delta1);...
           (delta3-omega^2*delta1);...
           (delta3-omega^2*delta1)].*(K*s1)-2*omega*[delta2;...
                                                  delta2;...
                                                   delta2].*(K*E*s1));

s2 = [sin(omega*t2);cos(omega*t2)];
                                               
                                           
pos2 = pc+ry*(K*s2);
acc2 = -ry*omega^2*K*s2;

delta1f = fliplr(delta1);
delta2f = fliplr(delta2);
delta3f = fliplr(delta3);

s3 = [sin(omega*t3);cos(omega*t3)];
pos3 = pc+ry*[delta1f;delta1f;delta1f].*(K*s3);
acc3 = ry*([(delta3f-omega^2*delta1f);...
           (delta3f-omega^2*delta1f);...
           (delta3f-omega^2*delta1f)].*(K*s3)-2*omega*[delta2f;...
                                                  delta2f;...
                                                   delta2f].*(K*E*s3));

pos = [pos1,pos2,pos3];
acc = [acc1,acc2,acc3];
t = [t1,t2,t3];
it1 = length(t1);
it2 = length([t1,t2]);
it3 = length(t);
% delta1y = ry*delta1;
% delta1x = rx*delta1;
% delta1z = rz*delta1;
% delta2x = rx*delta2;
% delta2y = ry*delta2;
% delta2z = rz*delta2;
% delta3x = delta3*rx;
% delta3y = delta3*ry;
% delta3z = delta3*rz;
% % setting the full trajectory
% delta1y = [delta1y,ry*ones(1,length(t2))];
% delta1x = [delta1x,rx*ones(1,length(t2))];
% delta1z = [delta1z,rz*ones(1,length(t2))];
% delta2x = [delta2x,zeros(1,length(t2))];
% delta2y = [delta2y,zeros(1,length(t2))];
% delta2z = [delta2z,zeros(1,length(t2))];
% delta3x = [delta3x,zeros(1,length(t2))];
% delta3y = [delta3y,zeros(1,length(t2))];
% delta3z = [delta3z,zeros(1,length(t2))];
% % setting the inverse transition trajectory
% delta1y = [delta1y,ry*fliplr(delta1)];
% delta1x = [delta1x,rx*fliplr(delta1)];
% delta1z = [delta1z,rz*fliplr(delta1)];
% delta2x = [delta2x,rx*fliplr(delta2)];
% delta2y = [delta2y,ry*fliplr(delta2)];
% delta2z = [delta2z,rz*fliplr(delta2)];
% delta3x = [delta3x,rx*fliplr(delta3)];
% delta3y = [delta3y,ry*fliplr(delta3)];
% delta3z = [delta3z,rz*fliplr(delta3)];


% vit(1:3,:)  = [delta2x.*sin(omega*t+px)+(omega.*delta1x.*cos(omega*t+px));
%                delta2y.*sin(omega*t+py)+(omega.*delta1y.*cos(omega*t+py));
%                delta2z.*sin(omega*t+pz)+(omega.*delta1z.*cos(omega*t+pz))];
% acc(1:3,:) = [(delta3x-delta1x.*omega^2).*sin(omega*t+px)+...
%     2*delta2x.*omega.*cos(omega*t+px);...
%     (delta3y-delta1y.*omega^2).*sin(omega*t+py)+...
%     2*delta2y.*omega.*cos(omega*t+py);
%     (delta3z-delta1z.*omega^2).*sin(omega*t+py)+...
%     2*delta2z.*omega.*cos(omega*t+pz)];
end

