function area = calc_area_SW(alpha,R,r,rc,pc)
%CALC_AREA_SW Calculates the area of the static workspace and returns it

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

pi.c(:,1) = [cx(1,1) cy(1,1);cx(1,2) cy(1,2)]\-[ck(1,1);ck(1,2)];
pi.c(:,2) = [cx(1,1) cy(1,1);cx(2,1) cy(2,1)]\-[ck(1,1);ck(2,1)];
pi.c(:,3) = [cx(1,1) cy(1,1);cx(2,2) cy(2,2)]\-[ck(1,1);ck(2,2)];
pi.c(:,4) = [cx(1,1) cy(1,1);cx(3,1) cy(3,1)]\-[ck(1,1);ck(3,1)];
pi.c(:,5) = [cx(1,1) cy(1,1);cx(3,2) cy(3,2)]\-[ck(1,1);ck(3,2)];
pi.c(:,6) = [cx(1,2) cy(1,2);cx(2,1) cy(2,1)]\-[ck(1,2);ck(2,1)];
pi.c(:,7) = [cx(1,2) cy(1,2);cx(2,2) cy(2,2)]\-[ck(1,2);ck(2,2)];
pi.c(:,8) = [cx(1,2) cy(1,2);cx(3,1) cy(3,1)]\-[ck(1,2);ck(3,1)];
pi.c(:,9) = [cx(1,2) cy(1,2);cx(3,2) cy(3,2)]\-[ck(1,2);ck(3,2)];
pi.c(:,10) = [cx(2,1) cy(2,1);cx(2,2) cy(2,2)]\-[ck(2,1);ck(2,2)];
pi.c(:,11) = [cx(2,1) cy(2,1);cx(3,1) cy(3,1)]\-[ck(2,1);ck(3,1)];
pi.c(:,12) = [cx(2,1) cy(2,1);cx(3,2) cy(3,2)]\-[ck(2,1);ck(3,2)];
pi.c(:,13) = [cx(2,2) cy(2,2);cx(3,1) cy(3,1)]\-[ck(2,2);ck(3,1)];
pi.c(:,14) = [cx(2,2) cy(2,2);cx(3,2) cy(3,2)]\-[ck(2,2);ck(3,2)];
pi.c(:,15) = [cx(3,1) cy(3,1);cx(3,2) cy(3,2)]\-[ck(3,1);ck(3,2)];


n =numel(fieldnames(pi));
possib = nchoosek(1:15,2);
% test = (1:length(possib(:,1)))';
possib = [(1:length(possib(:,1)))',possib];

t = linspace(0,1,10);
% Calcul des points milieu de chaque point d'intersection
lim_segs = [];
for i =1:length(possib(:,1))
    [seg.mp(i,:),seg.len(i),seg.v(i,:)] = segment_calc(pi.c(:,possib(i,2)),pi.c(:,possib(i,3)));
    seg.p1(:,i) = pi.c(:,possib(i,2));
    seg.p2(:,i) = pi.c(:,possib(i,3));
    return_bool(i) = check_cond(seg.p1(:,i),seg.p2(:,i),R,r,rc,alpha,pc);
    if return_bool(i) == 1
        lim_segs = [lim_segs,i];
    end
end

xseg = []; yseg = [];
for i =1:length(lim_segs)
    xseg =[xseg,pi.c(1,possib(lim_segs(i),2)),pi.c(1,possib(lim_segs(i),3))];
    yseg =[yseg,pi.c(2,possib(lim_segs(i),2)),pi.c(2,possib(lim_segs(i),3))];
end
seg1 = [xseg;yseg];
seg2 = unique(seg1','rows');
if ~isempty(seg2)
    [seg_X,seg_Y] = vec_acw_order(seg2(:,1),seg2(:,2));
    seg = [seg_X, seg_Y];
    area = polyarea(seg(:,1),seg(:,2));
else
    area = 0;
end

end

