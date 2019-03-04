function [fig1_h] = Robot_setup(R,r,alpha,rc,phic,hc)
%ROBOT_SETUP prints two images of the robot in the setup given by the
%inputs

% Plotting the first figure
fig1_h = figure;
hold on;
circle_plot(R,[0,0,0],'-b');








    function circle_h = circle_plot(radius,center,tnc)
        res = 200;
        if size(center) == [1,3]
            center = center';
        end
        t= linspace(0,2*pi,res);
        circle = [center(1);center(2)] + [radius*cos(t);...
                                          radius*sin(t);...
                                          ones(1,length(t))*center(3)];
                                     
        circle_h = plot3(circle(1,:),circle(2,:),circle(3,:),tnc);
    end
end

