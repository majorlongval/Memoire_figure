function [fig1_h,fig2_h] = Robot_setup(R,r,alpha,rc,phic,hc)
%ROBOT_SETUP prints two images of the robot in the setup given by the
%inputs

% Plotting the first figure

p = [0; 0; 2];
a11 = r*[0; 1; 0];
a12 = -a11;
ct = cos(2*pi/3); st = sin(2*pi/3);
Qs = [ct -st 0;st ct 0;0 0 1];
a21 = Qs*a11;
a22 = -a21;
a31 = Qs*a21;
a32 = -a31;
R1 = R*[cos(alpha);sin(alpha);0];
R2 = Qs*R1;
R3 = Qs*R2;

t = linspace(0,1,10); % parameter for length of cable

% Defining the lines
line11 = R1+a11 + (p-R1).*t;
line12 = R1+a12 + (p-R1).*t;
line21 = R2+a21 + (p-R2).*t;
line22 = R2+a22 + (p-R2).*t;
line31 = R3+a31 + (p-R3).*t;
line32 = R3+a32 + (p-R3).*t;

figure;
plot3(line11(1,:),line11(2,:),line11(3,:),'-k');
hold on;grid on;
plot3(line21(1,:),line21(2,:),line21(3,:),'--k');
plot3(line31(1,:),line31(2,:),line31(3,:),'-.k');
circle_plot(R,[0,0,0],'-b');  
circle_plot(r,[0,0,2],'-r');
plot3(line12(1,:),line12(2,:),line12(3,:),'-k');
plot3([line11(1,1),line12(1,1)],[line11(2,1),line12(2,1)],[line11(3,1),line12(3,1)],'-k');
plot3(line22(1,:),line22(2,:),line22(3,:),'--k');
plot3([line21(1,1),line22(1,1)],[line21(2,1),line22(2,1)],[line21(3,1),line22(3,1)],'--k');
plot3(line32(1,:),line32(2,:),line32(3,:),'-.k');
plot3([line31(1,1),line32(1,1)],[line31(2,1),line32(2,1)],[line31(3,1),line32(3,1)],'-.k');
plot3(rc*cos(phic),rc*sin(phic),2+hc,'.g','MarkerSize',15);
legend('legend1','legend2','legend3');
view(0,90);
fig1_h = gca;

figure;
plot3(line11(1,:),line11(2,:),line11(3,:),'-k');
grid on;
hold on;
plot3(line21(1,:),line21(2,:),line21(3,:),'--k');
plot3(line31(1,:),line31(2,:),line31(3,:),'-.k');
circle_plot(R,[0,0,0],'-b');  
circle_plot(r,[0,0,2],'-r');
plot3(line12(1,:),line12(2,:),line12(3,:),'-k');
plot3(line22(1,:),line22(2,:),line22(3,:),'--k');
plot3(line32(1,:),line32(2,:),line32(3,:),'-.k');
plot3(rc*cos(phic),rc*sin(phic),2+hc,'.g','MarkerSize',15);
view(0,0);
legend('legend1','legend2','legend3');
fig2_h = gca;







    function circle_h = circle_plot(radius,center,tnc)
        res = 200;
        if size(center) == [1,3]
            center = center';
        end
        t= linspace(0,2*pi,res);
        circle = [center(1);center(2);0] + [radius*cos(t);...
                                          radius*sin(t);...
                                          ones(1,length(t))*center(3)];
                                     
        circle_h = plot3(circle(1,:),circle(2,:),circle(3,:),tnc);
    end
end

