function [K,vol,l_seg_new] = calc_vol_INTER(alpha, R, r, r_c, phi_c, ...
    h_c, tx, ty, tz, m, Mx,...
    My, Mz,max_z,discret_z)
% CALC_VOL_INTER Calculates the volume of the intersection between the 
% static workspace when an external wrench is applied and the static 
% workspace when no wrench is applied 
%    alpha: construction angle
%       R : big radius
%       r : small radius
%       rc: distance of CM from centre
%       pc: angle of CM from x axis
%       hc: height of CM rel to centre
%       tx: x component of wrench
%       ty: y component of wrench
%       tz: z component of wrench
% Mx,My,Mz: Moment components
%       m : total mass
%    max_z: The maximum value of z evaluated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gravity
g = 9.81;
% Wrench component in the direction of gravity
Tz = m*g + tz;

% The discretization z
zr = linspace(0,max_z,discret_z);
dz = zr(2)-zr(1);


%% Calculating the coefficients of each line of the ETST and limits
[ay,az,ak,bx,bz,bk,cx,cy,ck,dx,dk,ey,ek,hz,hk] = ...
    calc_coeff(alpha, R, r, r_c, phi_c, h_c);
A = zeros(length(zr),6);
B = A;
C = A;
for j =1:length(zr)
    for i =1:6
        A(i,j) = ty*bx(i) + Tz*cx(i) + Mx*dx(i);
        B(i,j) = tx*ay(i) + Tz*cy(i) + My*ey(i);
        C(i,j) = tx*(az(i)*zr(j)+ak(i))+ty*(bz(i)*zr(j)+bk(i))+Tz*ck(i)+...
            Mx*dk(i)+My*ek(i)+Mz*(hz(i)*zr(j)+hk(i));
    end
end

%% Calculating all the intersections 
vol = 0;
l_seg_new = [];
for j = 1:(length(zr))
    inter(1,:,1,j) = [A(1,j) B(1,j);A(2,j) B(2,j)]\-[C(1,j);C(2,j)];
    inter(2,:,1,j) = [A(1,j) B(1,j);A(3,j) B(3,j)]\-[C(1,j);C(3,j)];
    inter(3,:,1,j) = [A(1,j) B(1,j);A(4,j) B(4,j)]\-[C(1,j);C(4,j)];
    inter(4,:,1,j) = [A(1,j) B(1,j);A(5,j) B(5,j)]\-[C(1,j);C(5,j)];
    inter(5,:,1,j) = [A(1,j) B(1,j);A(6,j) B(6,j)]\-[C(1,j);C(6,j)];
    inter(6,:,1,j) = [A(1,j) B(1,j);cx(1)  cy(1) ]\-[C(1,j);ck(1) ];
    inter(7,:,1,j) = [A(1,j) B(1,j);cx(2)  cy(2) ]\-[C(1,j);ck(2) ];
    inter(8,:,1,j) = [A(1,j) B(1,j);cx(3)  cy(3) ]\-[C(1,j);ck(3) ];
    inter(9,:,1,j) = [A(1,j) B(1,j);cx(4)  cy(4) ]\-[C(1,j);ck(4) ];
    inter(10,:,1,j) = [A(1,j) B(1,j);cx(5)  cy(5) ]\-[C(1,j);ck(5) ];
    inter(11,:,1,j) = [A(1,j) B(1,j);cx(6)  cy(6) ]\-[C(1,j);ck(6) ];
    
    inter(1,:,2,j) = [A(2,j) B(2,j);A(1,j) B(1,j)]\-[C(2,j);C(1,j)];
    inter(2,:,2,j) = [A(2,j) B(2,j);A(3,j) B(3,j)]\-[C(2,j);C(3,j)];
    inter(3,:,2,j) = [A(2,j) B(2,j);A(4,j) B(4,j)]\-[C(2,j);C(4,j)];
    inter(4,:,2,j) = [A(2,j) B(2,j);A(5,j) B(5,j)]\-[C(2,j);C(5,j)];
    inter(5,:,2,j) = [A(2,j) B(2,j);A(6,j) B(6,j)]\-[C(2,j);C(6,j)];
    inter(6,:,2,j) = [A(2,j) B(2,j);cx(1)  cy(1) ]\-[C(2,j);ck(1) ];
    inter(7,:,2,j) = [A(2,j) B(2,j);cx(2)  cy(2) ]\-[C(2,j);ck(2) ];
    inter(8,:,2,j) = [A(2,j) B(2,j);cx(3)  cy(3) ]\-[C(2,j);ck(3) ];
    inter(9,:,2,j) = [A(2,j) B(2,j);cx(4)  cy(4) ]\-[C(2,j);ck(4) ];
    inter(10,:,2,j) = [A(2,j) B(2,j);cx(5)  cy(5) ]\-[C(2,j);ck(5) ];
    inter(11,:,2,j) = [A(2,j) B(2,j);cx(6)  cy(6) ]\-[C(2,j);ck(6) ];
    
    inter(1,:,3,j) = [A(3,j) B(3,j);A(1,j) B(1,j)]\-[C(3,j);C(1,j)];
    inter(2,:,3,j) = [A(3,j) B(3,j);A(2,j) B(2,j)]\-[C(3,j);C(2,j)];
    inter(3,:,3,j) = [A(3,j) B(3,j);A(4,j) B(4,j)]\-[C(3,j);C(4,j)];
    inter(4,:,3,j) = [A(3,j) B(3,j);A(5,j) B(5,j)]\-[C(3,j);C(5,j)];
    inter(5,:,3,j) = [A(3,j) B(3,j);A(6,j) B(6,j)]\-[C(3,j);C(6,j)];
    inter(6,:,3,j) = [A(3,j) B(3,j);cx(1)  cy(1) ]\-[C(3,j);ck(1) ];
    inter(7,:,3,j) = [A(3,j) B(3,j);cx(2)  cy(2) ]\-[C(3,j);ck(2) ];
    inter(8,:,3,j) = [A(3,j) B(3,j);cx(3)  cy(3) ]\-[C(3,j);ck(3) ];
    inter(9,:,3,j) = [A(3,j) B(3,j);cx(4)  cy(4) ]\-[C(3,j);ck(4) ];
    inter(10,:,3,j) = [A(3,j) B(3,j);cx(5)  cy(5) ]\-[C(3,j);ck(5) ];
    inter(11,:,3,j) = [A(3,j) B(3,j);cx(6)  cy(6) ]\-[C(3,j);ck(6) ];
    
    inter(1,:,4,j) = [A(4,j) B(4,j);A(1,j) B(1,j)]\-[C(4,j);C(1,j)];
    inter(2,:,4,j) = [A(4,j) B(4,j);A(2,j) B(2,j)]\-[C(4,j);C(2,j)];
    inter(3,:,4,j) = [A(4,j) B(4,j);A(3,j) B(3,j)]\-[C(4,j);C(3,j)];
    inter(4,:,4,j) = [A(4,j) B(4,j);A(5,j) B(5,j)]\-[C(4,j);C(5,j)];
    inter(5,:,4,j) = [A(4,j) B(4,j);A(6,j) B(6,j)]\-[C(4,j);C(6,j)];
    inter(6,:,4,j) = [A(4,j) B(4,j);cx(1)  cy(1) ]\-[C(4,j);ck(1) ];
    inter(7,:,4,j) = [A(4,j) B(4,j);cx(2)  cy(2) ]\-[C(4,j);ck(2) ];
    inter(8,:,4,j) = [A(4,j) B(4,j);cx(3)  cy(3) ]\-[C(4,j);ck(3) ];
    inter(9,:,4,j) = [A(4,j) B(4,j);cx(4)  cy(4) ]\-[C(4,j);ck(4) ];
    inter(10,:,4,j) = [A(4,j) B(4,j);cx(5)  cy(5) ]\-[C(4,j);ck(5) ];
    inter(11,:,4,j) = [A(4,j) B(4,j);cx(6)  cy(6) ]\-[C(4,j);ck(6) ];
    
    inter(1,:,5,j) = [A(5,j) B(5,j);A(1,j) B(1,j)]\-[C(5,j);C(1,j)];
    inter(2,:,5,j) = [A(5,j) B(5,j);A(2,j) B(2,j)]\-[C(5,j);C(2,j)];
    inter(3,:,5,j) = [A(5,j) B(5,j);A(3,j) B(3,j)]\-[C(5,j);C(3,j)];
    inter(4,:,5,j) = [A(5,j) B(5,j);A(4,j) B(4,j)]\-[C(5,j);C(4,j)];
    inter(5,:,5,j) = [A(5,j) B(5,j);A(6,j) B(6,j)]\-[C(5,j);C(6,j)];
    inter(6,:,5,j) = [A(5,j) B(5,j);cx(1)  cy(1) ]\-[C(5,j);ck(1) ];
    inter(7,:,5,j) = [A(5,j) B(5,j);cx(2)  cy(2) ]\-[C(5,j);ck(2) ];
    inter(8,:,5,j) = [A(5,j) B(5,j);cx(3)  cy(3) ]\-[C(5,j);ck(3) ];
    inter(9,:,5,j) = [A(5,j) B(5,j);cx(4)  cy(4) ]\-[C(5,j);ck(4) ];
    inter(10,:,5,j) = [A(5,j) B(5,j);cx(5)  cy(5) ]\-[C(5,j);ck(5) ];
    inter(11,:,5,j) = [A(5,j) B(5,j);cx(6)  cy(6) ]\-[C(5,j);ck(6) ];
    
    inter(1,:,6,j) = [A(6,j) B(6,j);A(1,j) B(1,j)]\-[C(6,j);C(1,j)];
    inter(2,:,6,j) = [A(6,j) B(6,j);A(2,j) B(2,j)]\-[C(6,j);C(2,j)];
    inter(3,:,6,j) = [A(6,j) B(6,j);A(3,j) B(3,j)]\-[C(6,j);C(3,j)];
    inter(4,:,6,j) = [A(6,j) B(6,j);A(4,j) B(4,j)]\-[C(6,j);C(4,j)];
    inter(5,:,6,j) = [A(6,j) B(6,j);A(5,j) B(5,j)]\-[C(6,j);C(5,j)];
    inter(6,:,6,j) = [A(6,j) B(6,j);cx(1)  cy(1) ]\-[C(6,j);ck(1) ];
    inter(7,:,6,j) = [A(6,j) B(6,j);cx(2)  cy(2) ]\-[C(6,j);ck(2) ];
    inter(8,:,6,j) = [A(6,j) B(6,j);cx(3)  cy(3) ]\-[C(6,j);ck(3) ];
    inter(9,:,6,j) = [A(6,j) B(6,j);cx(4)  cy(4) ]\-[C(6,j);ck(4) ];
    inter(10,:,6,j) = [A(6,j) B(6,j);cx(5)  cy(5) ]\-[C(6,j);ck(5) ];
    inter(11,:,6,j) = [A(6,j) B(6,j);cx(6)  cy(6) ]\-[C(6,j);ck(6) ];
    
    inter(1,:,7,j)  =  [cx(1)  cy(1);A(1,j) B(1,j)]\-[ck(1);C(1,j)];
    inter(2,:,7,j)  =  [cx(1)  cy(1);A(2,j) B(2,j)]\-[ck(1);C(2,j)];
    inter(3,:,7,j)  =  [cx(1)  cy(1);A(3,j) B(3,j)]\-[ck(1);C(3,j)];
    inter(4,:,7,j)  =  [cx(1)  cy(1);A(4,j) B(4,j)]\-[ck(1);C(4,j)];
    inter(5,:,7,j)  =  [cx(1)  cy(1);A(5,j) B(5,j)]\-[ck(1);C(5,j)];
    inter(6,:,7,j)  =  [cx(1)  cy(1);A(6,j) B(6,j)]\-[ck(1);C(6,j)];
    inter(7,:,7,j)  =  [cx(1)  cy(1);cx(2)  cy(2) ]\-[ck(1);ck(2) ];
    inter(8,:,7,j)  =  [cx(1)  cy(1);cx(3)  cy(3) ]\-[ck(1);ck(3) ];
    inter(9,:,7,j)  =  [cx(1)  cy(1);cx(4)  cy(4) ]\-[ck(1);ck(4) ];
    inter(10,:,7,j) =  [cx(1)  cy(1);cx(5)  cy(5) ]\-[ck(1);ck(5) ];
    inter(11,:,7,j) =  [cx(1)  cy(1);cx(6)  cy(6) ]\-[ck(1);ck(6) ];

    inter(1,:,8,j)  =  [cx(2)  cy(2);A(1,j) B(1,j)]\-[ck(2);C(1,j)];
    inter(2,:,8,j)  =  [cx(2)  cy(2);A(2,j) B(2,j)]\-[ck(2);C(2,j)];
    inter(3,:,8,j)  =  [cx(2)  cy(2);A(3,j) B(3,j)]\-[ck(2);C(3,j)];
    inter(4,:,8,j)  =  [cx(2)  cy(2);A(4,j) B(4,j)]\-[ck(2);C(4,j)];
    inter(5,:,8,j)  =  [cx(2)  cy(2);A(5,j) B(5,j)]\-[ck(2);C(5,j)];
    inter(6,:,8,j)  =  [cx(2)  cy(2);A(6,j) B(6,j)]\-[ck(2);C(6,j)];
    inter(7,:,8,j)  =  [cx(2)  cy(2);cx(1)  cy(1) ]\-[ck(2);ck(1) ];
    inter(8,:,8,j)  =  [cx(2)  cy(2);cx(3)  cy(3) ]\-[ck(2);ck(3) ];
    inter(9,:,8,j)  =  [cx(2)  cy(2);cx(4)  cy(4) ]\-[ck(2);ck(4) ];
    inter(10,:,8,j) =  [cx(2)  cy(2);cx(5)  cy(5) ]\-[ck(2);ck(5) ];
    inter(11,:,8,j) =  [cx(2)  cy(2);cx(6)  cy(6) ]\-[ck(2);ck(6) ];
    
    inter(1,:,9,j)  =  [cx(3)  cy(3);A(1,j) B(1,j)]\-[ck(3);C(1,j)];
    inter(2,:,9,j)  =  [cx(3)  cy(3);A(2,j) B(2,j)]\-[ck(3);C(2,j)];
    inter(3,:,9,j)  =  [cx(3)  cy(3);A(3,j) B(3,j)]\-[ck(3);C(3,j)];
    inter(4,:,9,j)  =  [cx(3)  cy(3);A(4,j) B(4,j)]\-[ck(3);C(4,j)];
    inter(5,:,9,j)  =  [cx(3)  cy(3);A(5,j) B(5,j)]\-[ck(3);C(5,j)];
    inter(6,:,9,j)  =  [cx(3)  cy(3);A(6,j) B(6,j)]\-[ck(3);C(6,j)];
    inter(7,:,9,j)  =  [cx(3)  cy(3);cx(1)  cy(1) ]\-[ck(3);ck(1) ];
    inter(8,:,9,j)  =  [cx(3)  cy(3);cx(2)  cy(2) ]\-[ck(3);ck(2) ];
    inter(9,:,9,j)  =  [cx(3)  cy(3);cx(4)  cy(4) ]\-[ck(3);ck(4) ];
    inter(10,:,9,j) =  [cx(3)  cy(3);cx(5)  cy(5) ]\-[ck(3);ck(5) ];
    inter(11,:,9,j) =  [cx(3)  cy(3);cx(6)  cy(6) ]\-[ck(3);ck(6) ];
    
    inter(1,:,10,j)  =  [cx(4)  cy(4);A(1,j) B(1,j)]\-[ck(4);C(1,j)];
    inter(2,:,10,j)  =  [cx(4)  cy(4);A(2,j) B(2,j)]\-[ck(4);C(2,j)];
    inter(3,:,10,j)  =  [cx(4)  cy(4);A(3,j) B(3,j)]\-[ck(4);C(3,j)];
    inter(4,:,10,j)  =  [cx(4)  cy(4);A(4,j) B(4,j)]\-[ck(4);C(4,j)];
    inter(5,:,10,j)  =  [cx(4)  cy(4);A(5,j) B(5,j)]\-[ck(4);C(5,j)];
    inter(6,:,10,j)  =  [cx(4)  cy(4);A(6,j) B(6,j)]\-[ck(4);C(6,j)];
    inter(7,:,10,j)  =  [cx(4)  cy(4);cx(1)  cy(1) ]\-[ck(4);ck(1) ];
    inter(8,:,10,j)  =  [cx(4)  cy(4);cx(2)  cy(2) ]\-[ck(4);ck(2) ];
    inter(9,:,10,j)  =  [cx(4)  cy(4);cx(3)  cy(3) ]\-[ck(4);ck(3) ];
    inter(10,:,10,j) =  [cx(4)  cy(4);cx(5)  cy(5) ]\-[ck(4);ck(5) ];
    inter(11,:,10,j) =  [cx(4)  cy(4);cx(6)  cy(6) ]\-[ck(4);ck(6) ];
    
    inter(1,:,11,j)  =  [cx(5)  cy(5);A(1,j) B(1,j)]\-[ck(5);C(1,j)];
    inter(2,:,11,j)  =  [cx(5)  cy(5);A(2,j) B(2,j)]\-[ck(5);C(2,j)];
    inter(3,:,11,j)  =  [cx(5)  cy(5);A(3,j) B(3,j)]\-[ck(5);C(3,j)];
    inter(4,:,11,j)  =  [cx(5)  cy(5);A(4,j) B(4,j)]\-[ck(5);C(4,j)];
    inter(5,:,11,j)  =  [cx(5)  cy(5);A(5,j) B(5,j)]\-[ck(5);C(5,j)];
    inter(6,:,11,j)  =  [cx(5)  cy(5);A(6,j) B(6,j)]\-[ck(5);C(6,j)];
    inter(7,:,11,j)  =  [cx(5)  cy(5);cx(1)  cy(2) ]\-[ck(5);ck(1) ];
    inter(8,:,11,j)  =  [cx(5)  cy(5);cx(2)  cy(3) ]\-[ck(5);ck(2) ];
    inter(9,:,11,j)  =  [cx(5)  cy(5);cx(3)  cy(3) ]\-[ck(5);ck(3) ];
    inter(10,:,11,j) =  [cx(5)  cy(5);cx(4)  cy(4) ]\-[ck(5);ck(4) ];
    inter(11,:,11,j) =  [cx(5)  cy(5);cx(6)  cy(6) ]\-[ck(5);ck(6) ];
    

    inter(1,:,12,j)  =  [cx(6)  cy(6);A(1,j) B(1,j)]\-[ck(6);C(1,j)];
    inter(2,:,12,j)  =  [cx(6)  cy(6);A(2,j) B(2,j)]\-[ck(6);C(2,j)];
    inter(3,:,12,j)  =  [cx(6)  cy(6);A(3,j) B(3,j)]\-[ck(6);C(3,j)];
    inter(4,:,12,j)  =  [cx(6)  cy(6);A(4,j) B(4,j)]\-[ck(6);C(4,j)];
    inter(5,:,12,j)  =  [cx(6)  cy(6);A(5,j) B(5,j)]\-[ck(6);C(5,j)];
    inter(6,:,12,j)  =  [cx(6)  cy(6);A(6,j) B(6,j)]\-[ck(6);C(6,j)];
    inter(7,:,12,j)  =  [cx(6)  cy(6);cx(1)  cy(2) ]\-[ck(6);ck(1) ];
    inter(8,:,12,j)  =  [cx(6)  cy(6);cx(2)  cy(3) ]\-[ck(6);ck(2) ];
    inter(9,:,12,j)  =  [cx(6)  cy(6);cx(3)  cy(3) ]\-[ck(6);ck(3) ];
    inter(10,:,12,j) =  [cx(6)  cy(6);cx(4)  cy(4) ]\-[ck(6);ck(4) ];
    inter(11,:,12,j) =  [cx(6)  cy(6);cx(5)  cy(5) ]\-[ck(6);ck(5) ];
    
    l_seg = [];
    l_seg_ord = [];
    for i =1:12
        [ord_inter(:,1,i,j),ord_inter(:,2,i,j)] = ...
            vec_acw_order_lin(inter(:,1,i,j),inter(:,2,i,j));
        for k =1:10
            bool = check_cond_INTER(ord_inter(k,:,i,j),ord_inter(k+1,:,i,j),...
                   alpha, R, r, r_c, phi_c, h_c,tx, ty, tz, m, Mx,My, Mz,zr(j));
        if bool == 1
            l_seg = [l_seg;ord_inter(k,:,i,j);ord_inter(k+1,:,i,j)];
        end
        end
    end
    l_seg = unique(l_seg,'rows');

    if ~isempty(l_seg)
    [l_seg_ord(:,1),l_seg_ord(:,2)] = vec_acw_order(l_seg(:,1),l_seg(:,2));
%     vol =vol + polyarea(l_seg_ord(:,1),l_seg_ord(:,2))*dz;
    l_seg_new = [l_seg_new;[l_seg_ord,ones(length(l_seg_ord(:,1)),1)*zr(j)]];
    %else
%        vol = vol + 0;
    end
%     l_seg_res{j} = l_seg_new;
     
end
[K,vol] = convhulln(l_seg_new);
end