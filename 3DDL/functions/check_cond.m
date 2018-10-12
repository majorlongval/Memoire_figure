function return_bool = check_cond(P1,P2,R,r,rc,alpha,pc)
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

t = linspace(0,1,10);
roundn = @(x,n) round(x.*10.^n)./10.^n;

for k =1:length(t)
    pos(:,k) = P1 + t(k)*(P2-P1);
    for j =1:2
        for i = 1:3
            cond(i,j,k) = roundn(cx(i,j)*pos(1,k)+cy(i,j)*pos(2,k)+ck(i,j),4)>=0;
        end
    end
end

if all(cond(:) == 1)
    for j =1:2
        for i = 1:3
            cond2(i,j) = roundn(cx(i,j)*pos(1,5)+cy(i,j)*pos(2,5)+ck(i,j),8)==0;
        end
    end
    if any(cond2(:) ==1)
        return_bool = 1;
    else
        return_bool = 0;
    end
else
    return_bool = 0;
end


