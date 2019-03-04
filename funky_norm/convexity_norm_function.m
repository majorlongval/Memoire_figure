%% Graph to show that the norm function if convexe 
clear all; close all; clc;

% Corners of the box

center = [5 5 5];
amp = linspace(0.1,2,5);
for i =1:length(amp)
    
amplitudes(:,:,i) = amp(i)*[1 -1 -1;
               -1 -1 -1;
               -1 -1 1;
                1 -1 1;
                1 1 1;
               -1 1 1;
               -1 1 -1;
                1 1 -1];
            corners = zeros(8,3);
    for j =1:length(amplitudes(:,1,i))
     corners(j,:) = center + amplitudes(j,:,i); 
    end
    lp = calc_lin_interpol(corners,100);
    for j =1:length(lp(:,1))
    lpn(j) = norm(lp(j,:));
    end
    plot(lpn);
    hold on;
end
         



function list_of_points = calc_lin_interpol(points,nb)
list_of_points = [];
for i =1:length(points(:,1))-1
   for j = 1:3 
       temp(:,j) = linspace(points(i,j),points(i+1,j),nb);
   end
   list_of_points = [list_of_points;temp];
end
end

