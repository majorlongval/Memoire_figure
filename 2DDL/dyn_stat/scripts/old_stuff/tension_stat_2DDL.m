%% Graphique des tensions dans l'espace statique
clear all; close all;clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
%1) Input des valeurs pour lignes standard
L = 1;
l = 0.2;
ay = -0.072;
cy = -0.144;
cx = 0;
g = 9.81;

%2) Input des valeurs pour lignes max
ay2 = 0;
cy2 = 0;

%3) Calcul des des Ai et Bi standard
A1 = -L;
A2 = ay*(cy-L)+2*L*cy-l*(ay+L);
A3 = ay*(L-cy)-2*L*cy-l*(ay+L);

B1 = 1;
B2 = ay-l;
B3 = -(ay+l);

%4) Calcul des des Ai et Bi max

A12 = -L;
A22 = ay2*(cy2-L)+2*L*cy2-l*(ay2+L);
A32 = ay2*(L-cy2)-2*L*cy2-l*(ay2+L);

B12 = 1;
B22 = ay2-l;
B32 = -(ay2+l);


% lignes standard
a1 = [0;ay];
a2 = [0;-l];
a3 = [0;l];

b1 = [0;-L];
b2 = [0;L-l];
b3 = [0;L+l];

c = [cx;cy];

% lignes max
a12 = [0;ay2];
a22 = [0;-l];
a32 = [0;l];

b12 = [0;-L];
b22 = [0;L-l];
b32 = [0;L+l];

c2 = [cx;cy2];


E = [0,-1;1,0];
% Vecteur des accélérations
gam = [g;0.5*g;0.6*g];

x = linspace(1,5,5);
y = linspace(-L,L,200);



for j = 1:length(x)
    for i =1:length(y)
        % standard
        e1 = ([x(j);y(i)]+a1-b1)/norm([x(j);y(i)]+a1-b1);
        e2 = ([x(j);y(i)]+a2-b2)/norm([x(j);y(i)]+a2-b2);
        delta1    = (a1-c)'*E*e1;
        delta2    = (a2-c)'*E*e2;
        delta3    = (a3-c)'*E*e2;
        epsi1     = cross([e2;delta2],[e2;delta3]);
        epsi2     = cross([e2;delta3],[e1;delta1]);
        epsi3     = cross([e1;delta1],[e2;delta2]);
        M         = [e1 e2 e2;delta1 delta2 delta3]; 
        tau(:,i,j) = M\gam;
        tau2(:,i,j) = [epsi1';epsi2';epsi3']*gam/det(M);
        tau2(:,i,j) = [epsi1';epsi2';epsi3']*gam/((a2-a3 )'*E*e2*e1'*E*e2);
        
        % max
        e12 = ([x(j);y(i)]+a12-b12)/norm([x(j);y(i)]+a12-b12);
        e22 = ([x(j);y(i)]+a22-b22)/norm([x(j);y(i)]+a22-b22);
        delta12    = (c2-a12)'*E*e12;
        delta22    = (c2-a22)'*E*e22;
        delta32    = (c2-a32)'*E*e22;
        epsi12     = cross([e22;delta22],[e22;delta32]);
        epsi22     = cross([e22;delta32],[e12;delta12]);
        epsi32     = cross([e12;delta12],[e22;delta22]);
        M2         = [e12 e22 e22;delta12 delta22 delta32]; 
        tau22(:,i,j) = [epsi12';epsi22';epsi32']*gam/det(M2);
        
        
    end
end

fig = figure
grid minor;
hold on;
%axis([-L/L,L/L,-1.5,1.5]);

for j =5
    %standard
    plot(y/L,tau2(1,:,j)/g,'-b');
    plot(y/L,tau2(2,:,j)/g,'-r');
    plot(y/L,tau2(3,:,j)/g,'-k');
    % max
%     plot(y/L,tau22(1,:,j)/g,'--b');
%     plot(y/L,tau22(3,:,j)/g,'--k');
end
% %standard
plot([-A3/B3/L,-A3/B3/L],[-15,15],'--c');
plot([-A1/B1/L,-A1/B1/L],[-15,15],'--m');
% plot([-A32/B32/L,-A32/B32/L],[-15,15],'--g');
% plot([-A2/B2/L,-A2/B2/L],[-15,15],'--r');
plot([-1,1],[0,0],':k');
% max

xticks(linspace(-L/L,L/L,11));

xlabel('s01'); 
ylabel('s02');
legend_h=legend({'f1','f2','f3','f12','f32'},'Location','southeast');

title('title');



text(-A3/B3/L,-4,'yeq');
print('/home/jordan/Documents/Maitrise/recherche/memoire_figures/2DDL/dyn_stat/figs/espace_stat_tension','-depsc2');

        
        





        