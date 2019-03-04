clear all; close all; clc;

x = 5;
y = linspace(-1,1,50);

L = 1;
l = linspace(0,2*L,100);

s1 =sqrt(2)*1./(2*l.*(2*L+l));
s2 = 3*l.^2+4*L^2;
s3 = 8*L^3-10*L.*l.^2;
s4 = 4*L^4-l.^2*L^2+l.^4;
for i =1:length(y)
    sigma1(:,i) = (s1/x).*(sqrt(s2.*(x^2+y(i)^2)+s3.*y(i)+s4));
    sigma2(:,i) = sqrt(2).*(1./(2*x.*l))*sqrt(x^2+(y(i)-L)^2);
end

plot(l,sigma2(:,5)-sigma1(:,5));
hold on;


% plot(l,sigma2(:,5));
% hold on;


figure
ymin = L*(5.*l.^2-4*L^2)./(4*L^2+3.*l.^2);
plot(l,ymin);