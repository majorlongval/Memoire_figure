function vol = calc_vol_ETST(alpha, R, r, r_c, phi_c, ...
    h_c, tx, ty, tz, m, Mx,...
    My, Mz,max_z)
%CALC_VOL_ETST Calculates the volume of the static workspace when an
%external wrench is applied on the end-effector
%   alpha: construction angle
%       R : big radius
%       r : small radius
%       rc: distance of CM from centre
%       pc: angle of CM from x axis
%       hc: height of CM rel to centre
%       tx: x component of wrench
%       ty: y component of wrench
%       tz: z component of wrench
% Mx,My,Mz: Moment components
%       m : mass total

g = 9.81;
Tz = m*g + tz;

zr = linspace(0,max_z,500);
dz = zr(2)-zr(1);
% Calculating the sins and cosines for faster calc speed

[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, r_c, phi_c, h_c);
for j =1:length(zr)
    for i =1:6
        A(i,j) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
        B(i,j) = tx*ay(i) + Tz*cy(i) + My*ey(i);
        C(i,j) = tx*(az(i)*zr(j)+ak(i))+ty*(bz(i)*zr(j)+bk(i))+Tz*ck(i)+...
            Mx*dk(i)+My*ek(i)+Mz*(hz(i)*zr(j)+hk(i));
    end
end
% Calcul de tous les points d'intersection
vol = 0;
for j = 1:(length(zr)-1)
    inter(1,:,1,j) = [A(1,j) B(1,j);A(2,j) B(2,j)]\-[C(1,j);C(2,j)];
    inter(2,:,1,j) = [A(1,j) B(1,j);A(3,j) B(3,j)]\-[C(1,j);C(3,j)];
    inter(3,:,1,j) = [A(1,j) B(1,j);A(4,j) B(4,j)]\-[C(1,j);C(4,j)];
    inter(4,:,1,j) = [A(1,j) B(1,j);A(5,j) B(5,j)]\-[C(1,j);C(5,j)];
    inter(5,:,1,j) = [A(1,j) B(1,j);A(6,j) B(6,j)]\-[C(1,j);C(6,j)];
    
    inter(1,:,2,j) = [A(2,j) B(2,j);A(1,j) B(1,j)]\-[C(2,j);C(1,j)];
    inter(2,:,2,j) = [A(2,j) B(2,j);A(3,j) B(3,j)]\-[C(2,j);C(3,j)];
    inter(3,:,2,j) = [A(2,j) B(2,j);A(4,j) B(4,j)]\-[C(2,j);C(4,j)];
    inter(4,:,2,j) = [A(2,j) B(2,j);A(5,j) B(5,j)]\-[C(2,j);C(5,j)];
    inter(5,:,2,j) = [A(2,j) B(2,j);A(6,j) B(6,j)]\-[C(2,j);C(6,j)];
    
    inter(1,:,3,j) = [A(3,j) B(3,j);A(1,j) B(1,j)]\-[C(3,j);C(1,j)];
    inter(2,:,3,j) = [A(3,j) B(3,j);A(2,j) B(2,j)]\-[C(3,j);C(2,j)];
    inter(3,:,3,j) = [A(3,j) B(3,j);A(4,j) B(4,j)]\-[C(3,j);C(4,j)];
    inter(4,:,3,j) = [A(3,j) B(3,j);A(5,j) B(5,j)]\-[C(3,j);C(5,j)];
    inter(5,:,3,j) = [A(3,j) B(3,j);A(6,j) B(6,j)]\-[C(3,j);C(6,j)];
    
    inter(1,:,4,j) = [A(4,j) B(4,j);A(1,j) B(1,j)]\-[C(4,j);C(1,j)];
    inter(2,:,4,j) = [A(4,j) B(4,j);A(2,j) B(2,j)]\-[C(4,j);C(2,j)];
    inter(3,:,4,j) = [A(4,j) B(4,j);A(3,j) B(3,j)]\-[C(4,j);C(3,j)];
    inter(4,:,4,j) = [A(4,j) B(4,j);A(5,j) B(5,j)]\-[C(4,j);C(5,j)];
    inter(5,:,4,j) = [A(4,j) B(4,j);A(6,j) B(6,j)]\-[C(4,j);C(6,j)];
    
    inter(1,:,5,j) = [A(5,j) B(5,j);A(1,j) B(1,j)]\-[C(5,j);C(1,j)];
    inter(2,:,5,j) = [A(5,j) B(5,j);A(2,j) B(2,j)]\-[C(5,j);C(2,j)];
    inter(3,:,5,j) = [A(5,j) B(5,j);A(3,j) B(3,j)]\-[C(5,j);C(3,j)];
    inter(4,:,5,j) = [A(5,j) B(5,j);A(4,j) B(4,j)]\-[C(5,j);C(4,j)];
    inter(5,:,5,j) = [A(5,j) B(5,j);A(6,j) B(6,j)]\-[C(5,j);C(6,j)];
    
    inter(1,:,6,j) = [A(6,j) B(6,j);A(1,j) B(1,j)]\-[C(6,j);C(1,j)];
    inter(2,:,6,j) = [A(6,j) B(6,j);A(2,j) B(2,j)]\-[C(6,j);C(2,j)];
    inter(3,:,6,j) = [A(6,j) B(6,j);A(3,j) B(3,j)]\-[C(6,j);C(3,j)];
    inter(4,:,6,j) = [A(6,j) B(6,j);A(4,j) B(4,j)]\-[C(6,j);C(4,j)];
    inter(5,:,6,j) = [A(6,j) B(6,j);A(5,j) B(5,j)]\-[C(6,j);C(5,j)];
    
    l_seg = [];
    l_seg_new = [];
    for i =1:6
        [ord_inter(:,1,i,j),ord_inter(:,2,i,j)] = ...
            vec_acw_order_lin(inter(:,1,i,j),inter(:,2,i,j));
        for k =1:4
            bool = check_COND_ETST(ord_inter(k,:,i,j),ord_inter(k+1,:,i,j),...
                   alpha, R, r, r_c, phi_c, h_c,tx, ty, tz, m, Mx,My, Mz,zr(j));
        if bool == 1
            l_seg = [l_seg;ord_inter(k,:,i,j);ord_inter(k+1,:,i,j)];
        end
        end
    end
    l_seg = unique(l_seg,'rows');
%
    if ~isempty(l_seg)
    [l_seg_new(:,1,j),l_seg_new(:,2,j)] = vec_acw_order(l_seg(:,1),l_seg(:,2));
%     plot3([l_seg_new(:,1,j);l_seg_new(1,1,j)],...
%           [l_seg_new(:,2,j);l_seg_new(1,2,j)],...
%           ones(length(l_seg_new(:,1,j))+1,1)*zr(j),'-r');
    hold on;
    vol =vol + polyarea(l_seg_new(:,1,j),l_seg_new(:,2,j))*dz;
    else
        vol = vol + 0;
    end
%     l_seg_res{j} = l_seg_new;
     
end
% set(gca,'Zdir','reverse');
end

