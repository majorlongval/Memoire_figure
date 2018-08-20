clear all;close all;clc;
addpath('functions');
l = 1;

L = 5;

pos = [linspace(-L,L,100);ones(1,100)];

ay = linspace(-l,l,100);
cy = linspace(-l,l,100);



for i =1:length(ay)
        B2(i)   = ay(i)-l; 
        B3(i)   = -(ay(i)+l);
    for j =1:length(cy)
        A3(i,j) = ay(i)*(L-cy(j))-2*L*cy(j)-l*(ay(i)+L);
        A2(i,j) = -ay(i)*(L-cy(j))+2*L*cy(j)-l*(ay(i)+L);
        S2(i,j) = (A2(i,j)/B2(i));
        if S2(i,j) ==inf || S2(i,j) ==-inf
            S2(i,j) = NaN;
        end
        
        S3(i,j) = (A3(i,j)/B3(i));
        if S3(i,j) ==inf || S3(i,j) ==-inf
            S3(i,j) = NaN;
        end
        sumtau = [];
        for k = 1:length(pos(1,:))
            tau(:,k) = tension_stat(pos(1,k),struct('L',L,'l',l,'ay',ay(i),...
                                                  'cy', cy(j), 'cx',0.1),'ax',0);
            sumtau = [sumtau,sumtaupos(tau(:,k))];
        end
        tau1 = tau(1,:);
        tau2 = tau(2,:);
        tau3 = tau(3,:);
        taum(i,j) = sum(tau2(tau2>=0))/k;
        taum2(i,j) = sum(sumtau)/k;
    end
end


surf(ay,cy,(L+min(S2,S3))/L);

figure 

surf(ay,cy,taum);

figure
surf(ay,cy,taum*(L+min(S2,S3))/L);
