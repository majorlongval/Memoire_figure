close all; clc;

DS = experiment_loader(EXPDYN2,1);




% Setting ythe experiment beginning
t_ini = 225.5;

% Setting the time to put the effector in position
t_up = t_ini%+27.5;%s

% Setting the intervals between each % lapse
t_lapse = 30; %s

t_50 = t_up;
t_60 = t_50+t_lapse;
t_70 = t_60+t_lapse;
t_80 = t_70+t_lapse;
t_90 = t_80+t_lapse;
t_100 = t_90+t_lapse;
% Setting the indices;
i_ini_50 = round((t_50/DS.time(end))*length(DS.time));
i_fin_50 = round((t_60/DS.time(end))*length(DS.time));

i_ini_60 = round((t_60/DS.time(end))*length(DS.time));
i_fin_60 = round((t_70/DS.time(end))*length(DS.time));

i_ini_70 = round((t_70/DS.time(end))*length(DS.time));
i_fin_70 = round((t_80/DS.time(end))*length(DS.time));

i_ini_80 = round((t_80/DS.time(end))*length(DS.time));
i_fin_80 = round((t_90/DS.time(end))*length(DS.time));

i_ini_90 = round((t_90/DS.time(end))*length(DS.time));
i_fin_90 = round((t_100/DS.time(end))*length(DS.time));

% figure;
% subplot
% plot(DS.time,DS.gyro_x);
% hold on;

% %% Powers gyro x
% 
%50

gyro_x_50 = fft(DS.gyro_x(i_ini_50:i_fin_50)-mean(DS.gyro_x(i_ini_50:i_fin_50)));
gyro_x_50 = fftshift(gyro_x_50);
n_50 = length(gyro_x_50);                    % number of samples
f_50 = (-n_50/2:(n_50/2)-1)*(DS.fs/n_50);    % frequency range
x_power_50 = abs(gyro_x_50).^2/n_50; 

%60
gyro_x_60 = fft(DS.gyro_x(i_ini_60:i_fin_60)-mean(DS.gyro_x(i_ini_60:i_fin_60)));
gyro_x_60 = fftshift(gyro_x_60);
n_60 = length(gyro_x_60);                    % number of samples
f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
x_power_60 = abs(gyro_x_60).^2/n_60; 



%70
gyro_x_70 = fft(DS.gyro_x(i_ini_70:i_fin_70)-mean(DS.gyro_x(i_ini_70:i_fin_70)));
gyro_x_70 = fftshift(gyro_x_70);
n_70 = length(gyro_x_70);                    % number of samples
f_70 = (-n_70/2:(n_70/2)-1)*(DS.fs/n_70);    % frequency range
x_power_70 = abs(gyro_x_70).^2/n_70;  

%80
gyro_x_80 = fft(DS.gyro_x(i_ini_80:i_fin_80)-mean(DS.gyro_x(i_ini_80:i_fin_80)));
gyro_x_80 = fftshift(gyro_x_80);
n_80 = length(gyro_x_80);                    % number of samples
f_80 = (-n_80/2:(n_80/2)-1)*(DS.fs/n_80);    % frequency range
x_power_80 = abs(gyro_x_80).^2/n_80; 

%90
gyro_x_90 = fft(DS.gyro_x(i_ini_90:i_fin_90)-mean(DS.gyro_x(i_ini_90:i_fin_90)));
gyro_x_90 = fftshift(gyro_x_90);
n_90 = length(gyro_x_90);                    % number of samples
f_90 = (-n_90/2:(n_90/2)-1)*(DS.fs/n_90);    % frequency range
x_power_90 = abs(gyro_x_90).^2/n_90; 



% %% Powers gyro y
% 
%50

gyro_y_50 = fft(DS.gyro_y(i_ini_50:i_fin_50)-mean(DS.gyro_y(i_ini_50:i_fin_50)));
gyro_y_50 = fftshift(gyro_y_50);
n_50 = length(gyro_y_50);                    % number of samples
f_50 = (-n_50/2:(n_50/2)-1)*(DS.fs/n_50);    % frequency range
y_power_50 = abs(gyro_y_50).^2/n_50; 

%60
gyro_y_60 = fft(DS.gyro_y(i_ini_60:i_fin_60)-mean(DS.gyro_y(i_ini_60:i_fin_60)));
gyro_y_60 = fftshift(gyro_y_60);
n_60 = length(gyro_y_60);                    % number of samples
f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
y_power_60 = abs(gyro_y_60).^2/n_60; 



%70
gyro_y_70 = fft(DS.gyro_y(i_ini_70:i_fin_70)-mean(DS.gyro_y(i_ini_70:i_fin_70)));
gyro_y_70 = fftshift(gyro_y_70);
n_70 = length(gyro_y_70);                    % number of samples
f_70 = (-n_70/2:(n_70/2)-1)*(DS.fs/n_70);    % frequency range
y_power_70 = abs(gyro_y_70).^2/n_70;  

%80
gyro_y_80 = fft(DS.gyro_y(i_ini_80:i_fin_80)-mean(DS.gyro_y(i_ini_80:i_fin_80)));
gyro_y_80 = fftshift(gyro_y_80);
n_80 = length(gyro_y_80);                    % number of samples
f_80 = (-n_80/2:(n_80/2)-1)*(DS.fs/n_80);    % frequency range
y_power_80 = abs(gyro_y_80).^2/n_80; 

%90
gyro_y_90 = fft(DS.gyro_y(i_ini_90:i_fin_90)-mean(DS.gyro_y(i_ini_90:i_fin_90)));
gyro_y_90 = fftshift(gyro_y_90);
n_90 = length(gyro_y_90);                    % number of samples
f_90 = (-n_90/2:(n_90/2)-1)*(DS.fs/n_90);    % frequency range
y_power_90 = abs(gyro_y_90).^2/n_90; 


% %% Powers gyro z
% 
%50

gyro_z_50 = fft(DS.gyro_z(i_ini_50:i_fin_50)-mean(DS.gyro_z(i_ini_50:i_fin_50)));
gyro_z_50 = fftshift(gyro_z_50);
n_50 = length(gyro_z_50);                    % number of samples
f_50 = (-n_50/2:(n_50/2)-1)*(DS.fs/n_50);    % frequency range
z_power_50 = abs(gyro_z_50).^2/n_50; 

%60
gyro_z_60 = fft(DS.gyro_z(i_ini_60:i_fin_60)-mean(DS.gyro_z(i_ini_60:i_fin_60)));
gyro_z_60 = fftshift(gyro_z_60);
n_60 = length(gyro_z_60);                    % number of samples
f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
z_power_60 = abs(gyro_z_60).^2/n_60; 



%70
gyro_z_70 = fft(DS.gyro_z(i_ini_70:i_fin_70)-mean(DS.gyro_z(i_ini_70:i_fin_70)));
gyro_z_70 = fftshift(gyro_z_70);
n_70 = length(gyro_z_70);                    % number of samples
f_70 = (-n_70/2:(n_70/2)-1)*(DS.fs/n_70);    % frequency range
z_power_70 = abs(gyro_z_70).^2/n_70;  

%80
gyro_z_80 = fft(DS.gyro_z(i_ini_80:i_fin_80)-mean(DS.gyro_z(i_ini_80:i_fin_80)));
gyro_z_80 = fftshift(gyro_z_80);
n_80 = length(gyro_z_80);                    % number of samples
f_80 = (-n_80/2:(n_80/2)-1)*(DS.fs/n_80);    % frequency range
z_power_80 = abs(gyro_z_80).^2/n_80; 

%90
gyro_z_90 = fft(DS.gyro_z(i_ini_90:i_fin_90)-mean(DS.gyro_z(i_ini_90:i_fin_90)));
gyro_z_90 = fftshift(gyro_z_90);
n_90 = length(gyro_z_90);                    % number of samples
f_90 = (-n_90/2:(n_90/2)-1)*(DS.fs/n_90);    % frequency range
z_power_90 = abs(gyro_z_90).^2/n_90; 




% Total freq power 
tot_power_50 = x_power_50+y_power_50+z_power_50;
tot_power_60 = x_power_60+y_power_60+z_power_60;
tot_power_70 = x_power_70+y_power_70+z_power_70;
tot_power_80 = x_power_80+y_power_80+z_power_80;
tot_power_90 = x_power_90+y_power_90+z_power_90;

% 
% 
% Plotting power;

figure
subplot(5,1,1)
plot(f_50,tot_power_50,'-b');grid on;
xlabel('freq'); ylabel('power');
axis([0 5 0 0.2])
subplot(5,1,2)
plot(f_60,tot_power_60,'-r');grid on;
xlabel('freq'); ylabel('power');
axis([0 5 0 0.2])
subplot(5,1,3)
plot(f_70,tot_power_70,'-k');grid on;
xlabel('freq'); ylabel('power');
axis([0 5 0 0.2])
subplot(5,1,4)
plot(f_80,tot_power_80,'-m');grid on;
xlabel('freq'); ylabel('power');
axis([0 5 0 0.2])
subplot(5,1,5)
plot(f_90,tot_power_90,'-g');grid on;
xlabel('freq'); ylabel('power');
axis([0 5 0 0.2])
% saveas(gca,'frequency_power_per','svg')
%Plotting the gyroscope

% 
% figure
% subplot(3,1,1);
% plot(DS.time(i_ini_50:i_fin_50),DS.gyro_x(i_ini_50:i_fin_50),'-b');
% grid on; hold on;
% plot(DS.time(i_ini_60:i_fin_60),DS.gyro_x(i_ini_60:i_fin_60),'-r');
% plot(DS.time(i_ini_70:i_fin_70),DS.gyro_x(i_ini_70:i_fin_70),'-k');
% plot(DS.time(i_ini_80:i_fin_80),DS.gyro_x(i_ini_80:i_fin_80),'-m');
% plot(DS.time(i_ini_90:i_fin_90),DS.gyro_x(i_ini_90:i_fin_90),'-g');
% xlabel('temps');
% ylabel('gyrox');
% hold on;
% axis([DS.time(i_ini_50) DS.time(i_fin_90) min(DS.gyro_x) max(DS.gyro_x)]);
% % for i =1:18
% %     plot([10*(i-1)+t_ini,10*(i-1)+t_ini],[-5,5],'-r');
% % end
% % axis([0 220 -5 5]);
% 
% 
% subplot(3,1,2)
% plot(DS.time(i_ini_50:i_fin_50),DS.gyro_y(i_ini_50:i_fin_50),'-b');
% grid on; hold on;
% plot(DS.time(i_ini_60:i_fin_60),DS.gyro_y(i_ini_60:i_fin_60),'-r');
% plot(DS.time(i_ini_70:i_fin_70),DS.gyro_y(i_ini_70:i_fin_70),'-k');
% plot(DS.time(i_ini_80:i_fin_80),DS.gyro_y(i_ini_80:i_fin_80),'-m');
% plot(DS.time(i_ini_90:i_fin_90),DS.gyro_y(i_ini_90:i_fin_90),'-g');
% xlabel('temps');
% ylabel('gyroy');
% hold on;
% axis([DS.time(i_ini_50) DS.time(i_fin_90) min(DS.gyro_y) max(DS.gyro_y)]);
% % for i =1:18
% %     plot([10*(i-1)+t_ini,10*(i-1)+t_ini],[-5,5],'-r');
% % end
% % axis([DS.time(i_ini_50) DS.time(i_fin_90) ...
% %      min(DS.gyro_y(i_ini_50:i_fin_50) 5]);
% 
% 
% subplot(3,1,3)
% plot(DS.time(i_ini_50:i_fin_50),DS.gyro_z(i_ini_50:i_fin_50),'-b');
% grid on; hold on;
% plot(DS.time(i_ini_60:i_fin_60),DS.gyro_z(i_ini_60:i_fin_60),'-r');
% plot(DS.time(i_ini_70:i_fin_70),DS.gyro_z(i_ini_70:i_fin_70),'-k');
% plot(DS.time(i_ini_80:i_fin_80),DS.gyro_z(i_ini_80:i_fin_80),'-m');
% plot(DS.time(i_ini_90:i_fin_90),DS.gyro_z(i_ini_90:i_fin_90),'-g');
% xlabel('temps');
% ylabel('gyroz');
% hold on;
% axis([DS.time(i_ini_50) DS.time(i_fin_90) min(DS.gyro_z) max(DS.gyro_z)]);
% % for i =1:18
% %     plot([10*(i-1)+t_ini,10*(i-1)+t_ini],[-5,5],'-r');
% % end
% % axis([0 220 -5 5]);

% saveas(gca,'gyro_xyz','svg');




% % Plotting the acceleration
% figure
% subplot(3,1,1)
% plot(DS.time,DS.accel_xf,'-b'); grid on;
% xlabel('temps');
% ylabel('accelx');
% hold on;
% for i =1:18
%     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% end
% 
% 
% subplot(3,1,2)
% plot(DS.time,DS.accel_yf,'-b'); grid on;
% xlabel('temps');
% ylabel('accely');
% hold on;
% for i =1:18
%     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% end
% 
% 
% subplot(3,1,3)
% plot(DS.time,DS.accel_zf,'-b'); grid on;
% xlabel('temps');
% ylabel('accelz');
% hold on;
% for i =1:18
%     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% end
% 
% 
% 
% % Plotting the orientation
% figure
% subplot(3,1,1)
% plot(DS.time,DS.roll,'-b');grid on;
% xlabel('temps');
% ylabel('roll');
% hold on;
% for i =1:18
%     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% end
% subplot(3,1,2)
% plot(DS.time,DS.pitch,'-b');grid on;
% xlabel('temps');
% ylabel('pitch');
% hold on;
% for i =1:18
%     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% end
% subplot(3,1,3)
% plot(DS.time,DS.yaw,'-b');grid on;
% xlabel('temps');
% ylabel('yaw');
% hold on;
% for i =1:18
%     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% end
% 
% 
% 
% 
