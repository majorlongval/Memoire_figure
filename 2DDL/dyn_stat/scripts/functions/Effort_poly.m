function [varargout] = Effort_poly(x, y, geo)
%Effort_poly calculates the polyhedron of force/torque that can be applied
%to the center of mass of the end-effector.
%   the first 7 inputs are positions and parameters
%   the last input tells if a plot should be generated
%   V is the volume of the poly
%   Mer is the range of the torque that can be applied at fy = fx = 0
%   fyr is the range of the latteral force that can be applied  at fx = Me
%   =0
%   (fxr is not calculated as it is always [-9.81, 0];
%   phandle is the handle of the plot
nOutputs = nargout;
varargout = cell(1,nOutputs);

g = 9.81;

[ay, cy, cx, L, l,A, B, C] = unpack_geo(geo);
% Derived values
A1 = -L;
A2 = ay*(cy-L)+2*L*cy-l*(ay+L);
A3 = ay*(L-cy)-2*L*cy-l*(ay+L);

B1 = 1;
B2 = ay-l;
B3 = -(ay+l);

D(1) = 0;
D(2) = 2*L+ay;
D(3) = -D(2);


% Calculating the vectors
u(:,1) = [[-(C(1)+B(1)*x) D(1); -(C(2)+B(2)*x) D(2)]\...
       (-g*[A(1)+B(1)*y;A(2)+B(2)*y]);0]-[0; 0; -g];
u(:,2) = [[-(C(1)+B(1)*x) D(1); -(C(3)+B(3)*x) D(3)]\...
       (-g*[A(1)+B(1)*y;A(3)+B(3)*y]);0]-[0; 0; -g];
u(:,3) = [[-(C(2)+B(2)*x) D(2); -(C(3)+B(3)*x) D(3)]\...
       (-g*[A(2)+B(2)*y;A(3)+B(3)*y]);0]-[0; 0; -g];


% Calculating the volume
varargout{1} = abs(det([u(:,1),u(:,2),u(:,3)]));

% Calculating the range of Me
varargout{2}= sort([-(A(2)+B(2)*y)*g/D(2),-(A(3)+B(3)*y)*g/D(3)]);

% Calculating the range of fy
varargout{3} = g*[min_abs(-(A(2)+B(2)*y)/(C(2)+B(2)*x),...
    -(A(3)+B(3)*y)/(C(3)+B(3)*x)),-(A(1)+B(1)*y)/(C(1)+B(1)*x)];

if nargout >=4
    c = [1 0 0;0 0 0; 1 0 0];
    
    varargout{4} = figure;
    plot3([0 u(1,1)]/g,[0 u(2,1)]/g,[-g 0]/g,'-r','LineWidth',2);
    hold on; grid on;
    patch([0 u(1,1) u(1,2)]/g,[0 u(2,1) u(2,2)]/g,[-g 0 0]/g,'green',...
        'FaceAlpha',0.3,'FaceVertexCData',c,'EdgeColor','flat',...
        'LineWidth',2);
    patch([0 u(1,1) u(1,3)]/g,[0 u(2,1) u(2,3)]/g,[-g 0 0]/g,'green',...
        'FaceAlpha',0.3,'FaceVertexCData',c,'EdgeColor','flat',...
        'LineWidth',2);
    patch([0 u(1,2) u(1,3)]/g,[0 u(2,2) u(2,3)]/g,[-g 0 0]/g,'green',...
        'FaceAlpha',0.3,'FaceVertexCData',c,'EdgeColor','flat',...
        'LineWidth',2);
    plot3([0 u(1,2)]/g,[0 u(2,2)]/g,[-g 0]/g,'-r','LineWidth',2);
    plot3([0 u(1,3)]/g,[0 u(2,3)]/g,[-g 0]/g,'-r','LineWidth',2);
    %plot3([u(1,1) u(1,2)]/g,[u(2,1) u(2,2)]/g,[0 0],'-k','LineWidth',2);
    %plot3([u(1,2) u(1,3)]/g,[u(2,2) u(2,3)]/g,[0 0],'-k','LineWidth',2);
    %plot3([u(1,1) u(1,3)]/g,[u(2,1) u(2,3)]/g,[0 0],'-k','LineWidth',2);
    plot3(u(1,1)/g,u(2,1)/g,0,'*b','LineWidth',2);
    plot3(u(1,2)/g,u(2,2)/g,0,'*b','LineWidth',2);
    plot3(u(1,3)/g,u(2,3)/g,0,'*b','LineWidth',2);
    
%    plot3([0 0],varargout{2},[0 0],'-b','LineWidth',2);
%    plot3(varargout{3},[0 0],[0 0],'-b','LineWidth',2);
    xlabel('001'); ylabel('002'); zlabel('003');
end
end

