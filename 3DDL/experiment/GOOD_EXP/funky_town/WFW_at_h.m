function pii = WFW_at_h(alpha, R, r, rc, phic, hc,m,wc,z)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Calculating all the possible planes
addpath('C:\Users\jorda\OneDrive\Desktop\Memoire\Memoire_figure\3DDL\functions');


% Defining the list of all projection lines at the plane z =zmin
L_l_i = [];
for i =1:length(wc(:,1))
    L_l_i = [L_l_i;Calc_lines(alpha, R, r, rc, phic, ...
        hc,m,z,wc(i,:))];
end

% Determining the points on the z=zmin plane which are good.
pii = [];
for i =1:length(L_l_i(:,1))-1
    for j=i+1:length(L_l_i(:,1))
        temp = Calc_intersection(L_l_i(i,:), L_l_i(j,:));
        bool = 1;
        for k =1:length(L_l_i(:,1))
            if k ~= i && k~= j
                verif = L_l_i(k,1)*temp(1)+L_l_i(k,2)*temp(2)+L_l_i(k,3) > 0;
                bool = bool*verif;
            end
        end
        if bool ==1
            pii = [pii;temp,z];
            
        end
    end
end
end

