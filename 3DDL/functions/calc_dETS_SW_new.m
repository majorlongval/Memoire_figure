function d_ETS = calc_dETS_SW_new(alpha,R,r,rc,pc)
% Calculating the sins and cosines for faster calc speed
sa = sin(alpha);
s2a = sin(2*alpha);
ca = cos(alpha);
c2a = cos(2*alpha);
spc = sin(pc);
cpc = cos(pc);

% Parameters of the linear functions
cx(1,1) = 2*r*ca^2+rc*spc;
cx(1,2) = 2*r*ca^2-rc*spc;
cx(2,1) = -r*ca^2+rc*spc-sqrt(3)*r*s2a/2;
cx(2,2) = -r*ca^2-rc*spc-sqrt(3)*r*s2a/2;
cx(3,1) = -r*ca^2+rc*spc+sqrt(3)*r*s2a/2;
cx(3,2) = -r*ca^2-rc*spc+sqrt(3)*r*s2a/2;

cy(1,1) = r*s2a-rc*cpc;
cy(1,2) = r*s2a+rc*cpc;
cy(2,1) = sqrt(3)*r*ca^2-r*s2a/2-rc*cpc;
cy(2,2) = sqrt(3)*r*ca^2-r*s2a/2+rc*cpc;
cy(3,1) = -sqrt(3)*r*ca^2-r*s2a/2-rc*cpc;
cy(3,2) = -sqrt(3)*r*ca^2-r*s2a/2+rc*cpc;

ck(1,1) = R*r*ca+2*R*rc*ca*spc;
ck(1,2) = R*r*ca-2*R*rc*ca*spc;
ck(2,1) = R*r*ca-R*rc*ca*spc-sqrt(3)*R*rc*ca*cpc;
ck(2,2) = R*r*ca+R*rc*ca*spc+sqrt(3)*R*rc*ca*cpc;
ck(3,1) = R*r*ca-R*rc*ca*spc+sqrt(3)*R*rc*ca*cpc;
ck(3,2) = R*r*ca+R*rc*ca*spc-sqrt(3)*R*rc*ca*cpc;


inter(:,1,1) = [cx(1,1) cy(1,1);cx(1,2) cy(1,2)]\-[ck(1,1);ck(1,2)];
inter(:,2,1) = [cx(1,1) cy(1,1);cx(2,1) cy(2,1)]\-[ck(1,1);ck(2,1)];
inter(:,3,1) = [cx(1,1) cy(1,1);cx(2,2) cy(2,2)]\-[ck(1,1);ck(2,2)];
inter(:,4,1) = [cx(1,1) cy(1,1);cx(3,1) cy(3,1)]\-[ck(1,1);ck(3,1)];
inter(:,5,1) = [cx(1,1) cy(1,1);cx(3,2) cy(3,2)]\-[ck(1,1);ck(3,2)];

inter(:,1,2) = [cx(1,2) cy(1,2);cx(1,1) cy(1,1)]\-[ck(1,2);ck(1,1)];
inter(:,2,2) = [cx(1,2) cy(1,2);cx(2,1) cy(2,1)]\-[ck(1,2);ck(2,1)];
inter(:,3,2) = [cx(1,2) cy(1,2);cx(2,2) cy(2,2)]\-[ck(1,2);ck(2,2)];
inter(:,4,2) = [cx(1,2) cy(1,2);cx(3,1) cy(3,1)]\-[ck(1,2);ck(3,1)];
inter(:,5,2) = [cx(1,2) cy(1,2);cx(3,2) cy(3,2)]\-[ck(1,2);ck(3,2)];

inter(:,1,3) = [cx(2,1) cy(2,1);cx(1,1) cy(1,1)]\-[ck(2,1);ck(1,1)];
inter(:,2,3)=  [cx(2,1) cy(2,1);cx(1,2) cy(1,2)]\-[ck(2,1);ck(1,2)];
inter(:,3,3) = [cx(2,1) cy(2,1);cx(2,2) cy(2,2)]\-[ck(2,1);ck(2,2)];
inter(:,4,3) = [cx(2,1) cy(2,1);cx(3,1) cy(3,1)]\-[ck(2,1);ck(3,1)];
inter(:,5,3) = [cx(2,1) cy(2,1);cx(3,2) cy(3,2)]\-[ck(2,1);ck(3,2)];

inter(:,1,4) = [cx(2,2) cy(2,2);cx(1,1) cy(1,1)]\-[ck(2,2);ck(1,1)];
inter(:,2,4) = [cx(2,2) cy(2,2);cx(1,2) cy(1,2)]\-[ck(2,2);ck(1,2)];
inter(:,3,4) = [cx(2,2) cy(2,2);cx(2,1) cy(2,1)]\-[ck(2,2);ck(2,1)];
inter(:,4,4) = [cx(2,2) cy(2,2);cx(3,1) cy(3,1)]\-[ck(2,2);ck(3,1)];
inter(:,5,4) = [cx(2,2) cy(2,2);cx(3,2) cy(3,2)]\-[ck(2,2);ck(3,2)];

inter(:,1,5) = [cx(3,1) cy(3,1);cx(1,1) cy(1,1)]\-[ck(3,1);ck(1,1)];
inter(:,2,5) = [cx(3,1) cy(3,1);cx(1,2) cy(1,2)]\-[ck(3,1);ck(1,2)];
inter(:,3,5) = [cx(3,1) cy(3,1);cx(2,1) cy(2,1)]\-[ck(3,1);ck(2,1)];
inter(:,4,5) = [cx(3,1) cy(3,1);cx(2,2) cy(2,2)]\-[ck(3,1);ck(2,2)];
inter(:,5,5) = [cx(3,1) cy(3,1);cx(3,2) cy(3,2)]\-[ck(3,1);ck(3,2)];

inter(:,1,6) = [cx(3,2) cy(3,2);cx(1,1) cy(1,1)]\-[ck(3,2);ck(1,1)];
inter(:,2,6) = [cx(3,2) cy(3,2);cx(1,2) cy(1,2)]\-[ck(3,2);ck(1,2)];
inter(:,3,6) = [cx(3,2) cy(3,2);cx(2,1) cy(2,1)]\-[ck(3,2);ck(2,1)];
inter(:,4,6) = [cx(3,2) cy(3,2);cx(2,2) cy(2,2)]\-[ck(3,2);ck(2,2)];
inter(:,5,6) = [cx(3,2) cy(3,2);cx(3,1) cy(3,1)]\-[ck(3,2);ck(3,1)];

for k =1:6
    [ord_inter(:,1,k),ord_inter(:,2,k)] = vec_acw_order_lin(inter(1,:,k)',...
                                                        inter(2,:,k)');
end
% col_ar = {'*r','*g','*b','*k','*m'};
% for i =1:length(ord_inter(:,1,1))
%     plot(ord_inter(i,1,1),ord_inter(i,2,1),col_ar{i});
%     hold on;
%     %plot(ord_inter(i,1,2),ord_inter(i,2,2),col_ar{i});
%     plot(ord_inter(i,1,3),ord_inter(i,2,3),col_ar{i});
% end
% hold on;
% plot(ord_inter(:,1,2),ord_inter(:,2,2),'*r');
% plot(ord_inter(:,1,3),ord_inter(:,2,3),'*r');
% plot(ord_inter(:,1,4),ord_inter(:,2,4),'*r');
% plot(ord_inter(:,1,5),ord_inter(:,2,5),'*r');
% plot(ord_inter(:,1,6),ord_inter(:,2,6),'*r');
l_seg = [];

for i = 1:6
     for j =1:4
     bool = check_cond_new(ord_inter(j,:,i),ord_inter(j+1,:,i),...
                                                         R,r,rc,alpha,pc);
        if bool == 1
            l_seg = [l_seg;ord_inter(j,:,i);ord_inter(j+1,:,i)];
        end
    end
end
l_seg = unique(l_seg,'rows');
%
if ~isempty(l_seg)
    [l_seg_new(:,1),l_seg_new(:,2)] = vec_acw_order(l_seg(:,1),l_seg(:,2));
    d_ETS_m = mean(l_seg_new);
    d_ETS = norm(d_ETS_m);
else
    d_ETS = 1000;
end



% %% Plotting stuff
% 
% xmax = 1.5;
% xmin = -1.5;
% 
% x = linspace(xmin,xmax,1000);
% 
% Qs = [cos(2*pi/3) -sin(2*pi/3) 0;
%     sin(2*pi/3)  cos(2*pi/3) 0;
%     0                     0  1];
% R1 = R*[cos(alpha);sin(alpha);0];
% R2 = Qs*R1;
% R3 = Qs*R2;
% 
% figure;
% hold on; grid on;
% axis([-1.2 1.2 -1.2 1.2]);
% axis square;
% set(gca,'Ydir','reverse')
% plot(R1(1),R1(2),'*r');
% plot(R2(1),R2(2),'*r');
% plot(R3(1),R3(2),'*r');
% circle(0,0,R);
% xlabel('x'); ylabel('y');
% 
% for i =1:3
%     for j =1:2
%         if cy(i,j) ~= 0
%             y(i,:,j) = (-cx(i,j)*x-ck(i,j))/cy(i,j);
%             h_y(i,j) = plot(x,y(i,:,j));
%         else
%         plot([-ck(i,j)/cx(i,j) -ck(i,j)/cx(i,j)],[-10 10]);
%         end
%     end
% end
% 
% h_y(1,1).Color = 'r'; h_y(1).LineStyle = '-';
% h_y(1,2).Color = 'r'; h_y(2).LineStyle = '--';
% h_y(2,1).Color = 'b'; h_y(3).LineStyle = '-';
% h_y(2,2).Color = 'b'; h_y(4).LineStyle = '--';
% h_y(3,1).Color = 'g'; h_y(5).LineStyle = '-';
% h_y(3,2).Color = 'g'; h_y(6).LineStyle = '--';
% grid on;
% 
% fh = fill([l_seg_new(:,1);l_seg_new(1,1)],...
%           [l_seg_new(:,2);l_seg_new(1,2)],[0.6,0.6,1]);
% fh.EdgeColor = 'none';
% fh.FaceAlpha = 0.5;











% function h = circle(x,y,r)
%     hold on
%     th = 0:pi/50:2*pi;
%     xunit = r * cos(th) + x;
%     yunit = r * sin(th) + y;
%     h = plot(xunit, yunit,'-k');
% end
end