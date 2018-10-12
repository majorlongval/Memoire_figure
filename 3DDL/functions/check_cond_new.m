function return_bool = check_cond_new(P1,P2,R,r,rc,alpha,pc)
%CHECK_COND verifies that all the conditions are respected for a point
s2a = sin(2*alpha);
ca = cos(alpha);
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
roundn = @(x,n) round(x.*10.^n)./10.^n;


for i = 1:3
    for j= 1:2
        if (roundn(cx(i,j)*P1(1)+cy(i,j)*P1(2)+ck(i,j),4)>=0) &&...
           (roundn(cx(i,j)*P2(1)+cy(i,j)*P2(2)+ck(i,j),4)>=0)
            list_elem(i,j) = 1;
        else
            list_elem(i,j) = 0;
        end
    end
end
if all(list_elem(:) == 1)
    return_bool = 1;
else
    return_bool = 0;
end