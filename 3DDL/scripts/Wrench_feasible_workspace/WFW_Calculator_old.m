function [vol,K,vertices] = WFW_Calculator(R,r,alpha,rc,pc,hc,m,wc,zr)
addpath('/home/jordan/Documents/Memoire_figure_git/3DDL/functions')
%WFW_Calculator Calculates the volume of the WFW on a given height
% m : mass
% wc : wc
% zr : range of z = [zmin,zmaz]; zmin different than zero.
zmin = zr(1);
zmax = zr(2);
phic = pc;
planes = [];
for k =1:length(wc(:,1))
    planes = [planes;Calc_planes(alpha,R,r,rc,phic,hc,m,...
        wc(k,:))];
end
% Adding extra planes for zmax and zmin
planes = [planes;0 0 -1,zmax;0 0 1 -zmin];



% Defining the list of all projection lines at the plane z =zmin
L_l_i = [];
for i =1:length(wc(:,1))
    L_l_i = [L_l_i;Calc_lines(alpha, R, r, rc, phic, ...
        hc,m,zmin,wc(i,:))];
end

% Determining the points on the z=zmin plane which are good.
p_inter_ini = [];

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
            p_inter_ini = [p_inter_ini;i,j,temp];
            
        end
    end
end

% From the list of points which was found previously, Determing if the
% planes which give rise to this point intersection with any other planes
n_p_o_l =[];
for i =1:length(p_inter_ini(:,1))
    for j =1:length(planes(:,1))-2
        if j~=p_inter_ini(i,1) && j~=p_inter_ini(i,2)
            temp = calc_plane_intersection(planes(p_inter_ini(i,1),:),...
                planes(p_inter_ini(i,2),:),...
                planes(j,:));
            if ~isempty(temp)
                bool =1;
                for h = 1:length(planes(:,1))
                    if h~=p_inter_ini(i,1) && h~=p_inter_ini(i,2) && h~=j
                        valid = (planes(h,1:3)*temp)+planes(h,4)>0;
                        bool = bool*valid;
                    end
                end
                if bool ==1
                    n_p_o_l =[n_p_o_l;p_inter_ini(i,1),p_inter_ini(i,2),h,temp'];
                end
            end
            
        end
    end
end

% Defining the list of all projection lines at the plane z =zmax
L_l_e = [];
for i =1:length(wc(:,1))
    L_l_e = [L_l_e;Calc_lines(alpha, R, r, rc, phic, ...
        hc,m,zmax,wc(i,:))];%,...
    %rhomax,rhomin,fmin,fmax)];
end

% Determining the points on the z=zmin plane which are good.
p_inter_end = [];
for i =1:length(L_l_e(:,1))-1
    for j=i+1:length(L_l_e(:,1))
        temp = Calc_intersection(L_l_e(i,:), L_l_e(j,:));
        bool = 1;
        for k =1:length(L_l_e(:,1))
            if k ~= i && k~= j
                verif = L_l_e(k,1)*temp(1)+L_l_e(k,2)*temp(2)+L_l_e(k,3) > 0;
                bool = bool*verif;
            end
        end
        if bool ==1
            p_inter_end = [p_inter_end;i,j,temp];
            
        end
    end
end

vertices = [p_inter_ini(:,3:4),ones(length(p_inter_ini(:,1)),1)*zmin;];

if ~isempty(n_p_o_l)
    vertices = [vertices;n_p_o_l(:,4:6)];
end

if ~isempty(p_inter_end)
    vertices = [vertices;...
        p_inter_end(:,3:4),ones(length(p_inter_end(:,1)),1)*zmax];
end
[K,vol] = convhulln(vertices);


    function [vert, set_new_planes] = branching(ip1,ip2,planes)
        vert = [];
        set_new_planes = [];
        for i =1:length(planes(:,1))-1
            if i ~= ip1 && i ~=ip2
                temp = calc_plane_intersection(planes(ip1,:),...
                    planes(ip2,:),...
                    planes(i,:));
                if ~isempty(temp)
                    bool =1;
                    for h = 1:length(planes(:,1))
                        if h~=ip1 && h~=ip2 && h~=i
                            valid = (planes(h,1:3)*temp)+planes(h,4)>0;
                            bool = bool*valid;
                        end
                    end
                    if bool ==1
                        vert = [vert;temp'];
                        set_new_planes = [set_new_planes;ip1,ip2,i];
                    end
                end
            end
        end
    end
end
