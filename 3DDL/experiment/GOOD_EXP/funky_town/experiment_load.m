close all; clc;
%Desired experimental data
WFW_at_h_script;
addpath('C:\Users\jorda\OneDrive\Desktop\Memoire\Memoire_figure\3DDL\functions');

set(groot, 'defaultFigureRenderer', 'painters')

% THE IMU DATA SHOULD BE LOADED. THE BEST WAY IS JUST TO DRAG AND DROP 
DS = experiment_loader(funkytownIMU,10);
% LOADING THE DATA FROM THE SIMULINK MODEL
load('funky_town.mat');


st_nf = 1;
st_nf_index = round(st_nf/opvar(1,end)*length((opvar(1,:))));

et_nf = 215;
et_nf_index = round((et_nf/opvar(1,end))*length((opvar(1,:))));



figure
plot3(opvar(2,st_nf_index:100:et_nf_index),-opvar(3,st_nf_index:100:et_nf_index),2*ones(1,length(opvar(4,st_nf_index:100:et_nf_index))),'-b');
K = convhull(pii_0(:,1),pii_0(:,2));
hold on; grid on;
plot3(pii_0(K,1),pii_0(K,2),pii_0(K,3),'-g');
%set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');

xlabel('x'); ylabel('y');

view(0,-90);

saveas(gca,'traj_mass','epsc');


st_nf = 55;
st_nf_index = round(st_nf/opvar(1,end)*length((opvar(1,:))));

figure
p1 = plot(opvar(1,st_nf_index:100:et_nf_index),opvar(2,st_nf_index:100:et_nf_index),'-b'); grid on; hold on;
p2 = plot(opvar(1,st_nf_index:100:et_nf_index),-opvar(3,st_nf_index:100:et_nf_index),'-r');
p3 = plot(opvar(1,st_nf_index:100:et_nf_index),opvar(4,st_nf_index:100:et_nf_index)+0.317,'-k');

st = 61;
for j =1:4
for i =1:5
    if i ~= 5
      duree = 3;
      sind = round(st/opvar(1,end)*length(opvar(1,:)));
      durind = round((duree)/opvar(1,end)*length(opvar(1,:)));
    if j == 4
      plot([opvar(1,sind),opvar(1,sind)],[-2 2.5],'--k');
      plot([opvar(1,sind+durind),opvar(1,sind+durind)],[-2 2.5],'--k'); 
    else
      plot([opvar(1,sind),opvar(1,sind)],[-2 2.5],'--k');
      plot([opvar(1,sind+durind),opvar(1,sind+durind)],[-2 2.5],'--k');
    end
    end
  st = st+8;
end
end

axis([55 215 -0.6 2.1]);
l1 = legend([p1,p2,p3],'x','y','z');
l1.Location = 'NorthEastOutside';
saveas(gca,'traj_mass_temps','epsc');

st = 10;
stindex = round((st/DS.time(end))*length(DS.time));

et = 170;
etindex = round((et/DS.time(end))*length(DS.time));


figure
p1 = plot(DS.time(stindex:etindex),DS.roll(stindex:etindex),'-b');hold on; grid on;
p2 = plot(DS.time(stindex:etindex),DS.pitch(stindex:etindex),'-r');
p3 = plot(DS.time(stindex:etindex),DS.yaw(stindex:etindex),'-k');

axis([10 170 -2 2]);

st = 14;
dv = {}
for j =1:4
for i =1:5
    if i ~= 5
[f,power,sind,durind] = fft_pallier(DS,st,3);
    dv{i,1,j} = f;
    dv{i,2,j} = power;
    dv{i,3,j} = sind;
    dv{i,4,j} = durind;
    if j == 4
      plot([DS.time(sind),DS.time(sind)],[-2 2],'--k');
  plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--k'); 
    else
  plot([DS.time(sind),DS.time(sind)],[-2 2],'--k');
  plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--k');
    end
    end
  st = st+8;
end
dv{i+1,1,j} = f;
temp  = [dv{1,2,j}';dv{2,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
dv{i+1,2,j} = mean(temp);
end
% axis([78 230 -2 1.5]);
xlabel('x'); ylabel('y');
l1 = legend([p1,p2,p3],'roll','pitch','yaw');
l1.Location = 'NorthEastOutside';

saveas(gca,'IMU_mass','epsc');

figure
for j = 1:4
    subplot(4,1,j)
    plot(dv{i+1,1,j},dv{i+1,2,j},'-b');
     axis([0 5 0 1]);
    grid on;
    xlabel('freq'); ylabel('mag');
    title(strcat('eta',num2str(j)));
end

saveas(gca,'IMU_mass_freq','epsc');






% Mass and force 


st = 220;
sti = round(st/opvar(1,end)*length((opvar(1,:))));
et = 380;
eti = round((et/opvar(1,end))*length((opvar(1,:))));

figure 
p1 = plot(opvar(1,sti:100:eti),opvar(2,sti:100:eti),'-b'); hold on; grid on;
p2 = plot(opvar(1,sti:100:eti),-opvar(3,sti:100:eti),'-r');
p3 = plot(opvar(1,sti:100:eti),opvar(4,sti:100:eti)+0.317,'-k');
axis([220 380 -0.5 2.1]);
xlabel('x'); ylabel('y');
st = 229;
for j =1:4
for i =1:5
    if i ~= 5
      duree = 3;
      sind = round(st/opvar(1,end)*length(opvar(1,:)));
      durind = round((duree)/opvar(1,end)*length(opvar(1,:)));
    if j == 4
      plot([opvar(1,sind),opvar(1,sind)],[-2 2.5],'--k');
      plot([opvar(1,sind+durind),opvar(1,sind+durind)],[-2 2.5],'--k'); 
    else
      plot([opvar(1,sind),opvar(1,sind)],[-2 2.5],'--k');
      plot([opvar(1,sind+durind),opvar(1,sind+durind)],[-2 2.5],'--k');
    end
    end
  st = st+8;
end
end

legend([p1,p2,p3],'x','y','z');

saveas(gca,'pos_time_mass_force','epsc');




figure
plot3(opvar(2,sti:100:eti),-opvar(3,sti:100:eti),opvar(4,sti:100:eti)+0.317,'-b');
grid on; hold on;
K = convhull(pii_1(:,1),pii_1(:,2));
plot3(pii_1(K,1),pii_1(K,2),2*ones(length(K)),'-g');
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
xlabel('x');
ylabel('y');
view(0,90);


saveas(gca,'pos_mass_force','epsc');


st = 185;
sti = round((st/DS.time(end))*length(DS.time));
et = 340;
eti = round((et/DS.time(end))*length(DS.time));

figure 
plot(DS.time(sti:eti),DS.roll(sti:eti),'-b');
hold on; grid on;
plot(DS.time(sti:eti),DS.pitch(sti:eti),'-r');
plot(DS.time(sti:eti),DS.yaw(sti:eti),'-k');
axis([st et -4 6]);
xlabel('x'); ylabel('y');

st = 190;
dv = {}
for j =1:4
for i =1:5
    if i ~= 5
[f,power,sind,durind] = fft_pallier(DS,st,3);
    dv{i,1,j} = f;
    dv{i,2,j} = power;
    dv{i,3,j} = sind;
    dv{i,4,j} = durind;
    if j == 4
      plot([DS.time(sind),DS.time(sind)],[-4 6],'--k');
  plot([DS.time(sind+durind),DS.time(sind+durind)],[-4 6],'--k'); 
    else
  plot([DS.time(sind),DS.time(sind)],[-4 6],'--k');
  plot([DS.time(sind+durind),DS.time(sind+durind)],[-4 6],'--k');
    end
    end
  st = st+8;
end
dv{i+1,1,j} = f;
if j ==1
temp =  [dv{1,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
else
temp  = [dv{1,2,j}';dv{2,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
end
dv{i+1,2,j} = mean(temp);
end
% axis([78 230 -2 1.5]);
xlabel('x'); ylabel('y');
l1 = legend([p1,p2,p3],'roll','pitch','yaw');
l1.Location = 'NorthEastOutside';

saveas(gca,'IMU_mass_force','epsc');

figure
for j = 1:4
    subplot(4,1,j)
    plot(dv{i+1,1,j},dv{i+1,2,j},'-b');
     axis([0 5 0 max(dv{i+1,2,j})]);
    grid on;
    xlabel('freq'); ylabel('mag');
    title(strcat('eta',num2str(j)));
end

saveas(gca,'IMU_mass_force_freq','epsc');



















































% 
% % NO FORCE
% 
% figure;
% plot(DS.time, DS.roll-mean(DS.roll)); hold on; grid on;
% plot(DS.time, DS.pitch-mean(DS.pitch));
% plot(DS.time, DS.yaw-mean(DS.yaw));
% % 
% st = 6;
% dv = {}
% for j =1:4
% for i =1:5
% [f,power,sind,durind] = fft_pallier(DS,st,4);
%     dv{i,1,j} = f;
%     dv{i,2,j} = power;
%     dv{i,3,j} = sind;
%     dv{i,4,j} = durind;
%     if j == 4
%       plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r'); 
%     else
%   plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r');
%     end
%   st = st+8;
% end
% dv{i+1,1,j} = f;
% temp  = [dv{1,2,j}';dv{2,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
% dv{i+1,2,j} = mean(temp);
% end
% figure
% for j = 1:4
%     subplot(5,1,j)
%     plot(dv{i+1,1,j},dv{i+1,2,j},'-b');
%     axis([0 5 0 max(dv{i+1,2,j})]);
% end
% 
% 
% figure;
% plot(DS.time, DS.roll-mean(DS.roll)); hold on; grid on;
% plot(DS.time, DS.pitch-mean(DS.pitch));
% plot(DS.time, DS.yaw-mean(DS.yaw));
% 
% st = 2;
% dvt = {}
% for j =1:4
% for i =1:5
% [f,power,sind,durind] = fft_pallier(DS,st,4);
%     dv{i,1,j} = f;
%     dv{i,2,j} = power;
%     dv{i,3,j} = sind;
%     dv{i,4,j} = durind;
%     if j == 4
%       plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r'); 
%     else
%   plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r');
%     end
%   st = st+8;
% end
% dv{i+1,1,j} = f;
% if j ==1
% temp  = [dv{1,2,j}';dv{2,2,j}';dv{4,2,j}';dv{5,2,j}'];
% else
%     temp  = [dv{1,2,j}';dv{2,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
% end
% dv{i+1,2,j} = mean(temp);
% end
% figure
% for j = 1:4
%     subplot(5,1,j)
%     plot(dv{i+1,1,j},dv{i+1,2,j},'-b');
%     axis([0 5 0 max(dv{i+1,2,j})]);
%     xlabel('frequency');
%     ylabel('power');
% end
% 
% 
% 
% 
% 
% % WITH FORCE
% figure;
% plot(DS.time, DS.roll-mean(DS.roll)); hold on; grid on;
% plot(DS.time, DS.pitch-mean(DS.pitch));
% plot(DS.time, DS.yaw-mean(DS.yaw));
% % 
% st = 182;
% dv = {}
% for j =1:4
% for i =1:5
% [f,power,sind,durind] = fft_pallier(DS,st,4);
%     dv{i,1,j} = f;
%     dv{i,2,j} = power;
%     dv{i,3,j} = sind;
%     dv{i,4,j} = durind;
%     if j == 4
%       plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r'); 
%     else
%   plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r');
%     end
%   st = st+8;
% end
% dv{i+1,1,j} = f;
% if j ==1
%     temp  = [dv{1,2,j}';dv{2,2,j}';dv{4,2,j}';dv{5,2,j}'];
% else
%     temp  = [dv{1,2,j}';dv{2,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
% end
% dv{i+1,2,j} = mean(temp);
% end
% figure
% for j = 1:4
%     subplot(5,1,j)
%     plot(dv{i+1,1,j},dv{i+1,2,j},'-b');
%     axis([0 5 0 max(dv{i+1,2,j})]);
% end
% % 
% % 
% figure;
% plot(DS.time, DS.roll-mean(DS.roll)); hold on; grid on;
% plot(DS.time, DS.pitch-mean(DS.pitch));
% plot(DS.time, DS.yaw-mean(DS.yaw));
% 
% st = 178;
% dvt = {}
% for j =1:4
% for i =1:5
% [f,power,sind,durind] = fft_pallier(DS,st,4);
%     dv{i,1,j} = f;
%     dv{i,2,j} = power;
%     dv{i,3,j} = sind;
%     dv{i,4,j} = durind;
%     if j == 4
%       plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r'); 
%     else
%   plot([DS.time(sind),DS.time(sind)],[-2 2],'--b');
%   plot([DS.time(sind+durind),DS.time(sind+durind)],[-2 2],'--r');
%     end
%   st = st+8;
% end
% dv{i+1,1,j} = f;
% if j ==1
%     temp  = [dv{1,2,j}';dv{2,2,j}';dv{4,2,j}';dv{5,2,j}'];
% else
%     temp  = [dv{1,2,j}';dv{2,2,j}';dv{3,2,j}';dv{4,2,j}';dv{5,2,j}'];
% end
% dv{i+1,2,j} = mean(temp);
% end
% figure
% for j = 1:4
%     subplot(5,1,j)
%     plot(dv{i+1,1,j},dv{i+1,2,j},'-b');
%     axis([0 5 0 max(dv{i+1,2,j})]);
%     xlabel('frequency');
%     ylabel('power');
% end










% figure
% plot(opvar(1,:),opvar(2,:)); grid on; hold on;
% plot(opvar(1,:),opvar(3,:));
% plot(opvar(1,:),opvar(4,:));
% 
% figure;
% plot(DS.time, DS.roll-mean(DS.roll)); hold on; grid on;
% plot(DS.time, DS.pitch-mean(DS.pitch));
% plot(DS.time, DS.yaw-mean(DS.yaw));
% 
% %% only_mass
% 
% [f_60_1,p_60_1] = fft_pallier(DS,3,8);
% [f_60_2,p_60_2] = fft_pallier(DS,3+8,8);
% [f_60_3,p_60_3] = fft_pallier(DS,3+16,8);
% [f_60_4,p_60_4] = fft_pallier(DS,3+24,8);
% [f_60_5,p_60_5] = fft_pallier(DS,3+32,8);
% P60 = [p_60_1';p_60_2';p_60_4';p_60_5'];
% P60m = mean(P60);
% 
% [f_70_1,p_70_1] = fft_pallier(DS,3+40,8);
% [f_70_2,p_70_2] = fft_pallier(DS,3+48,8);
% [f_70_3,p_70_3] = fft_pallier(DS,3+56,8);
% [f_70_4,p_70_4] = fft_pallier(DS,3+64,8);
% [f_70_5,p_70_5] = fft_pallier(DS,3+72,8);
% P70 = [p_70_1';p_70_2';p_70_3';p_70_4';p_70_5'];
% P70m = mean(P70);
% % 
% [f_80_1,p_80_1] = fft_pallier(DS,3+80,8);
% [f_80_2,p_80_2] = fft_pallier(DS,3+88,8);
% [f_80_3,p_80_3] = fft_pallier(DS,3+96,8);
% [f_80_4,p_80_4] = fft_pallier(DS,3+104,8);
% [f_80_5,p_80_5] = fft_pallier(DS,3+112,8);
% P80 = [p_80_1';p_80_2';p_80_3';p_80_4';p_80_5'];
% P80m = mean(P80);
% 
% [f_90_1,p_90_1] = fft_pallier(DS,3+120,8);
% [f_90_2,p_90_2] = fft_pallier(DS,3+128,8);
% [f_90_3,p_90_3] = fft_pallier(DS,3+136,8);
% [f_90_4,p_90_4] = fft_pallier(DS,3+144,8);
% [f_90_5,p_90_5] = fft_pallier(DS,3+152,8);
% P90 = [p_90_1';p_90_2';p_90_3';p_90_4';p_90_5'];
% P90m = mean(P90);

% %% Adding force
% 
% [f_60_1f,p_60_1f] = fft_pallier(DS,3+160,8);
% [f_60_2f,p_60_2f] = fft_pallier(DS,3+168,8);
% [f_60_3f,p_60_3f] = fft_pallier(DS,3+16,8);
% [f_60_4f,p_60_4f] = fft_pallier(DS,3+24,8);
% [f_60_5f,p_60_5f] = fft_pallier(DS,3+32,8);
% P60 = [p_60_1';p_60_2';p_60_3';p_60_4';p_60_5'];
% P60m = mean(P60);
% 
% [f_70_1,p_70_1] = fft_pallier(DS,3+40,8);
% [f_70_2,p_70_2] = fft_pallier(DS,3+48,8);
% [f_70_3,p_70_3] = fft_pallier(DS,3+56,8);
% [f_70_4,p_70_4] = fft_pallier(DS,3+64,8);
% [f_70_5,p_70_5] = fft_pallier(DS,3+72,8);
% P70 = [p_70_1';p_70_2';p_70_3';p_70_4';p_70_5'];
% P70m = mean(P70);
% % 
% [f_80_1,p_80_1] = fft_pallier(DS,3+80,8);
% [f_80_2,p_80_2] = fft_pallier(DS,3+88,8);
% [f_80_3,p_80_3] = fft_pallier(DS,3+96,8);
% [f_80_4,p_80_4] = fft_pallier(DS,3+104,8);
% [f_80_5,p_80_5] = fft_pallier(DS,3+112,8);
% P80 = [p_80_1';p_80_2';p_80_3';p_80_4';p_80_5'];
% P80m = mean(P80);
% 
% [f_90_1,p_90_1] = fft_pallier(DS,3+120,8);
% [f_90_2,p_90_2] = fft_pallier(DS,3+128,8);
% [f_90_3,p_90_3] = fft_pallier(DS,3+136,8);
% [f_90_4,p_90_4] = fft_pallier(DS,3+144,8);
% [f_90_5,p_90_5] = fft_pallier(DS,3+152,8);
% P90 = [p_90_1';p_90_2';p_90_3';p_90_4';p_90_5'];
% P90m = mean(P90);

% [f_60_1,p_60_1] = fft_pallier(DS,3,8);
% [f_60_2,p_60_2] = fft_pallier(DS,3+8,8);
% [f_60_3,p_60_3] = fft_pallier(DS,3+16,8);
% [f_60_4,p_60_4] = fft_pallier(DS,3+24,8);
% [f_60_5,p_60_5] = fft_pallier(DS,3+32,8);
% P60 = [p_60_1';p_60_2';p_60_3';p_60_4';p_60_5'];
% P60m = mean(P60);
% % 
% indices_just_mass = plot_IMU_roll_pitch_yaw(DS,20,8,5,4);
% fft_plot(indices_just_mass,DS);
% 
% indices_whole_chabang = plot_IMU_roll_pitch_yaw(DS,180,8,5,4);
% fft_plot(indices_whole_chabang,DS);
% 











function indices = plot_IMU_roll_pitch_yaw(DS,stime,inttime,nbpoints,nbcycles)
%% This function gives a nice plot of the data and splits the data
tDATA = DS.time;
rDATA = DS.roll;
pDATA = DS.pitch;
yDATA = DS.yaw;
    indices = stime*round(length(tDATA)/tDATA(end));
    for i =1:nbcycles
        indices = [indices,(stime+(i*nbpoints*inttime))*round(length(tDATA)/tDATA(end))];
    end
    min_data = min([rDATA(indices(1):indices(end));...
                                        pDATA(indices(1):indices(end));...
                                        yDATA(indices(1):indices(end))]);
    max_data = max([rDATA(indices(1):indices(end));...
                                        pDATA(indices(1):indices(end));...
                                        yDATA(indices(1):indices(end))]);
    figure;
    plot(tDATA,rDATA,'-b',tDATA,pDATA,'-r',tDATA,yDATA,'-k');
    hold on;
    for i =1:nbcycles
        plot([tDATA(indices(i)),tDATA(indices(i))],[min_data,max_data],'--k');
    end
     grid on;
    xlabel('time [s]');
    ylabel('angle [deg]');
    l1 = legend('roll','pitch','yaw');
    l1.Location  = 'best';
    axis([stime tDATA(indices(end)) min_data,...
                                    max_data]);
    
end

function fft_plot(indices,DS)
    figure;
    for i =1:length(indices)-1
        subplot(4,1,i);
        rd = DS.roll(indices(i):indices(i+1));
        rp = DS.pitch(indices(i):indices(i+1));
        ry = DS.yaw(indices(i):indices(i+1));
        norm = sqrt(rd.^2+rp.^2+ry.^2);
        norm = detrend(norm,'linear');
        fs = DS.fs;
        fftv = fft(norm);
        n = length(norm);
        f = (0:n-1)*(fs/n);
        power = abs(fftv).^2/n;
        plot(f,power);
        xlabel('Frequency');
        ylabel('Power')
        axis([0,1,0,max(power)]);
    end
% 
end

function [f,power,sind,durind] = fft_pallier(DS,stime,duree)
 sind = round(stime/DS.time(end)*length(DS.time));
        durind = round((duree)/DS.time(end)*length(DS.time));
        rd = DS.roll(sind:sind+durind)-DS.roll(sind);
        rp = DS.pitch(sind:sind+durind)-DS.pitch(sind);
        ry = DS.yaw(sind:sind+durind)-DS.yaw(sind);
        
       rd = detrend(rd,'linear'); rp = detrend(rp,'linear'); ry = detrend(ry,'linear');
        fftrd = fft(rd); fftrp = fft(rp); fftry = fft(ry);
%         norm = sqrt(rd.^2+rp.^2+ry.^2);
%         norm = norm-norm(1);
%         norm = norm-mean(norm);
%         norm = detrend(norm,'linear');
        fs = DS.fs;
%         fftv = fft(norm);
        n = length(rd);
        f = (0:n-1)*(fs/n);
        powerrd = abs(fftrd).^2/n%.^2/n;
        powerrp = abs(fftrp).^2/n%.^2/n;
        powerry = abs(fftry).^2/n%.^2/n;
        power = powerrd+powerrp+powerry;

end



% %% Only mass
% % Getting the indexes for only mass
% i_60d = round((20/DS.time(end))*length(DS.time));
% i_60e = round((57.5/DS.time(end))*length(DS.time));
% i_70d = i_60e;
% i_70e = round((95/DS.time(end))*length(DS.time));
% i_80d = i_70e;
% i_80e = round((132.5/DS.time(end))*length(DS.time));
% i_90d = i_80e;
% i_90e = round((170/DS.time(end))*length(DS.time));
% 
% % Plotting the graph of only mass
% figure;
% plot(DS.time(i_60d:i_90e),DS.roll(i_60d:i_90e),'-b'); % roll
% grid on; hold on;
% plot(DS.time(i_60d:i_90e),DS.pitch(i_60d:i_90e),'-r'); % pitch
% plot(DS.time(i_60d:i_90e),DS.yaw(i_60d:i_90e)-mean(DS.yaw(i_60d:i_90e)),'-k'); % yaw
% minval = -4; maxval = 4;
% axis([DS.time(i_60d) DS.time(i_90e) minval maxval]);
% plot([DS.time(i_60d),DS.time(i_60d)],[minval,maxval],'--k');
% plot([DS.time(i_70d),DS.time(i_70d)],[minval,maxval],'--k');
% plot([DS.time(i_80d),DS.time(i_80d)],[minval,maxval],'--k');
% plot([DS.time(i_90d),DS.time(i_90d)],[minval,maxval],'--k');
% 
% % Calculating the norm of the theta vector
% norm_60 = sqrt((DS.roll(i_60d:i_60e).^2)+...
%                (DS.pitch(i_60d:i_60e).^2)+...
%                (DS.yaw(i_60d:i_60e).^2));
% norm_60 = detrend(norm_60,'linear');
% % 
% norm_70 = sqrt((DS.roll(i_70d:i_70e).^2)+...
%                (DS.pitch(i_70d:i_70e).^2)+...
%                (DS.yaw(i_70d:i_70e).^2));
% norm_70 = detrend(norm_70,'linear');
% 
% norm_80 = sqrt((DS.roll(i_80d:i_80e).^2)+...
%                (DS.pitch(i_80d:i_80e).^2)+...
%                (DS.yaw(i_80d:i_80e).^2));
% norm_80 = detrend(norm_80,'linear');
% 
% norm_90 = sqrt((DS.roll(i_90d:i_90e).^2)+...
%                (DS.pitch(i_90d:i_90e).^2)+...
%                (DS.yaw(i_90d:i_90e).^2));
% norm_90 = detrend(norm_90,'linear');
%            
% % Doing the spectral analysis
% fs = 100;
% fft_norm_60 = fft(norm_60);
% n60 = length(norm_60);          % number of samples
% f60 = (0:n60-1)*(fs/n60);     % frequency range
% power60 = abs(fft_norm_60).^2/n60;    % power of the DFT
% 
% fft_norm_70 = fft(norm_70);
% n70 = length(norm_70);          % number of samples
% f70 = (0:n70-1)*(fs/n70);     % frequency range
% power70 = abs(fft_norm_70).^2/n70;    % power of the DFT
% 
% fft_norm_80 = fft(norm_80);
% n80 = length(norm_80);          % number of samples
% f80 = (0:n80-1)*(fs/n80);     % frequency range
% power80 = abs(fft_norm_80).^2/n80;    % power of the DFT
% 
% fft_norm_90 = fft(norm_90);
% n90 = length(norm_90);          % number of samples
% f90 = (0:n90-1)*(fs/n90);     % frequency range
% power90 = abs(fft_norm_90).^2/n90;    % power of the DFT
% 
% % Plotting what the
% 
% %
% % 
% % % 
% % % fft_roll_70 = fft(DS.roll(i_70d:i_70e));
% % % n70 = length(DS.roll(i_70d:i_70e));          % number of samples
% % % f70 = (0:n70-1)*(fs/n70);     % frequency range
% % % power70 = abs(fft_roll_70).^2/n70;    % power of the DFT
% % % 
% % % fft_roll_80 = fft(DS.roll(i_80d:i_80e));
% % % n80 = length(DS.roll(i_80d:i_80e));          % number of samples
% % % f80 = (0:n80-1)*(fs/n80);     % frequency range
% % % power80 = abs(fft_roll_80).^2/n80;    % power of the DFT
% % % 
% % % fft_roll_90 = fft(DS.roll(i_90d:i_90e));
% % % n90 = length(DS.roll(i_90d:i_90e));          % number of samples
% % % f90 = (0:n90-1)*(fs/n90);     % frequency range
% % % power90 = abs(fft_roll_90).^2/n90;    % power of the DFT
% % 
% figure;
% subplot(4,1,1)
% plot(f60(1:end),power60(1:end))
% xlabel('Frequency')
% ylabel('Power')
% axis([0,1,0,500]);
% subplot(4,1,2)
% plot(f70(1:end),power70(1:end))
% xlabel('Frequency')
% ylabel('Power')
% axis([0,1,0,500]);
% subplot(4,1,3)
% plot(f80(1:end),power80(1:end))
% xlabel('Frequency')
% ylabel('Power')
% axis([0,1,0,500]);
% subplot(4,1,4)
% plot(f90(1:end),power90(1:end))
% xlabel('Frequency')
% ylabel('Power')
% axis([0,1,0,500]);
% % subplot(4,1,2)
% % plot(f70,power70)
% % xlabel('Frequency')
% % ylabel('Power')
% % subplot(4,1,3);
% % plot(f80,power80)
% % xlabel('Frequency')
% % ylabel('Power')
% % subplot(4,1,4);
% % plot(f90,power90)
% % xlabel('Frequency')
% % ylabel('Power')
% 
% 
% % FINDING THE STARTING TIME
% k  = find(opvar(2,:)~=0);
% i_deb = k(1);
% t_deb = opvar(1,i_deb);
% 






%% Getting the data for only mass


% % RETRIEVING THE VALUES of the indices %
% 
% i_bary = i_deb + round(7.5*length(opvar(1,:))/opvar(1,end));
% i_60d = i_bary + round(7.5*length(opvar(1,:))/opvar(1,end));
% t_60d = opvar(1,i_60d);
% i_60e = i_60d + round(30*length(opvar(1,:))/opvar(1,end));
% t_60e = opvar(1,i_60e);
% i_70d = i_60e + round(7.5*length(opvar(1,:))/opvar(1,end));
% t_70d = opvar(1,i_70d);
% i_70e = i_70d + round(30*length(opvar(1,:))/opvar(1,end));
% t_70e = opvar(1,i_70e);
% i_80d = i_70e + round(7.5*length(opvar(1,:))/opvar(1,end));
% t_80d = opvar(1,i_80d);
% i_80e = i_80d + round(30*length(opvar(1,:))/opvar(1,end));
% t_80e = opvar(1,i_80e);
% i_90d = i_80e + round(7.5*length(opvar(1,:))/opvar(1,end));
% t_90d = opvar(1,i_90d);
% i_90e = i_90d + round(30*length(opvar(1,:))/opvar(1,end));
% t_90e = opvar(1,i_90e);
% 
% figure;
% grid on; hold on;
% % Plotting centre and barycentre
% plot3(opvar(2,1),opvar(3,1),opvar(4,1),'ok');
% plot3([opvar(2,1),opvar(2,i_bary)],[opvar(3,1),opvar(3,i_bary)],[opvar(4,1),opvar(4,i_bary)],'--k');
% plot3(opvar(2,i_bary),opvar(3,i_bary),opvar(4,i_bary),'ob');
% plot3([opvar(2,i_bary),opvar(2,i_60d)],[opvar(3,i_bary),opvar(3,i_60d)],[opvar(4,i_bary),opvar(4,i_60d)],'--k');
% % Plotting 60 %
% plot3(opvar(2,i_60d:i_60e),opvar(3,i_60d:i_60e),opvar(4,i_60d:i_60e),'-b');
% for i =2:5
% plot3(points_0(i,1),-points_0(i,2),points(i,3),'*b');
% end
% plot3([opvar(2,i_60e),opvar(2,i_70d)],[opvar(3,i_60e),opvar(3,i_70d)],[opvar(4,i_60e),opvar(4,i_70d)],'--k');
% % Plotting 70 %
% plot3(opvar(2,i_70d:i_70e),opvar(3,i_70d:i_70e),opvar(4,i_70d:i_70e),'-r');
% for i =7:11
% plot3(points_0(i,1),-points_0(i,2),points(i,3),'*r');
% end
% plot3([opvar(2,i_70e),opvar(2,i_80d)],[opvar(3,i_70e),opvar(3,i_80d)],[opvar(4,i_70e),opvar(4,i_80d)],'--k');
% % Plotting 80 %
% plot3(opvar(2,i_80d:i_80e),opvar(3,i_80d:i_80e),opvar(4,i_80d:i_80e),'-g');
% for i =12:16
% plot3(points_0(i,1),-points_0(i,2),points(i,3),'*g');
% end
% plot3([opvar(2,i_80e),opvar(2,i_90d)],[opvar(3,i_80e),opvar(3,i_90d)],[opvar(4,i_80e),opvar(4,i_90d)],'--k');
% % % % Plotting graphs from OPVAR
% % % 
% % figure; grid on;
% % plot3(opvar(2,:),opvar(3,:),opvar(4,:));
% % 
% % 
% % 
% figure;
% subplot(2,1,1)
% hold on;
% grid on;
% plot(opvar(1,:),opvar(2,:));
% plot(opvar(1,:),opvar(3,:));
% plot(opvar(1,:),opvar(4,:));

% 
% 
% 
% DS = experiment_loader(funkytownIMU,1);
% 
% subplot(2,1,2);
% plot(DS.time,DS.roll);
% hold on; grid on;
% plot(DS.time,DS.pitch);
% plot(DS.time,DS.yaw);
% for i =1:5
%     plot([10+40*(i-1),10+40*(i-1)],[-2,2],'--k');
% end
% 
% i_60_deb = round(10*length(DS.time)/DS.time(end));
% i_70_deb = round((50*length(DS.time)/DS.time(end)));
% 
% roll_60 = fft(DS.roll(i_60_deb:i_70_deb));
% n_60 = length(roll_60);                    % number of samples
% f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
% roll_60_power = abs(roll_60).^2/n_60; 
% f = (0:n-1)*(fs/n);     % frequency range
% power = abs(y).^2/n;    % power of the DFT
% % 
% % Setting ythe experiment beginning
% t_ini = 225.5;
% 
% % Setting the time to put the effector in position
% t_up = t_ini%+27.5;%s
% 
% % Setting the intervals between each % lapse
% t_lapse = 30; %s
% 
% t_50 = t_up;
% t_60 = t_50+t_lapse;
% t_70 = t_60+t_lapse;
% t_80 = t_70+t_lapse;
% t_90 = t_80+t_lapse;
% t_100 = t_90+t_lapse;
% % Setting the indices;
% i_ini_50 = round((t_50/DS.time(end))*length(DS.time));
% i_fin_50 = round((t_60/DS.time(end))*length(DS.time));
% 
% i_ini_60 = round((t_60/DS.time(end))*length(DS.time));
% i_fin_60 = round((t_70/DS.time(end))*length(DS.time));
% 
% i_ini_70 = round((t_70/DS.time(end))*length(DS.time));
% i_fin_70 = round((t_80/DS.time(end))*length(DS.time));
% 
% i_ini_80 = round((t_80/DS.time(end))*length(DS.time));
% i_fin_80 = round((t_90/DS.time(end))*length(DS.time));
% 
% i_ini_90 = round((t_90/DS.time(end))*length(DS.time));
% i_fin_90 = round((t_100/DS.time(end))*length(DS.time));
% 
% % figure;
% % subplot
% % plot(DS.time,DS.gyro_x);
% % hold on;
% 
% % %% Powers gyro x
% % 
% %50
% 
% gyro_x_50 = fft(DS.gyro_x(i_ini_50:i_fin_50)-mean(DS.gyro_x(i_ini_50:i_fin_50)));
% gyro_x_50 = fftshift(gyro_x_50);
% n_50 = length(gyro_x_50);                    % number of samples
% f_50 = (-n_50/2:(n_50/2)-1)*(DS.fs/n_50);    % frequency range
% x_power_50 = abs(gyro_x_50).^2/n_50; 
% 
% %60
% gyro_x_60 = fft(DS.gyro_x(i_ini_60:i_fin_60)-mean(DS.gyro_x(i_ini_60:i_fin_60)));
% gyro_x_60 = fftshift(gyro_x_60);
% n_60 = length(gyro_x_60);                    % number of samples
% f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
% x_power_60 = abs(gyro_x_60).^2/n_60; 
% 
% 
% 
% %70
% gyro_x_70 = fft(DS.gyro_x(i_ini_70:i_fin_70)-mean(DS.gyro_x(i_ini_70:i_fin_70)));
% gyro_x_70 = fftshift(gyro_x_70);
% n_70 = length(gyro_x_70);                    % number of samples
% f_70 = (-n_70/2:(n_70/2)-1)*(DS.fs/n_70);    % frequency range
% x_power_70 = abs(gyro_x_70).^2/n_70;  
% 
% %80
% gyro_x_80 = fft(DS.gyro_x(i_ini_80:i_fin_80)-mean(DS.gyro_x(i_ini_80:i_fin_80)));
% gyro_x_80 = fftshift(gyro_x_80);
% n_80 = length(gyro_x_80);                    % number of samples
% f_80 = (-n_80/2:(n_80/2)-1)*(DS.fs/n_80);    % frequency range
% x_power_80 = abs(gyro_x_80).^2/n_80; 
% 
% %90
% gyro_x_90 = fft(DS.gyro_x(i_ini_90:i_fin_90)-mean(DS.gyro_x(i_ini_90:i_fin_90)));
% gyro_x_90 = fftshift(gyro_x_90);
% n_90 = length(gyro_x_90);                    % number of samples
% f_90 = (-n_90/2:(n_90/2)-1)*(DS.fs/n_90);    % frequency range
% x_power_90 = abs(gyro_x_90).^2/n_90; 
% 
% 
% 
% % %% Powers gyro y
% % 
% %50
% 
% gyro_y_50 = fft(DS.gyro_y(i_ini_50:i_fin_50)-mean(DS.gyro_y(i_ini_50:i_fin_50)));
% gyro_y_50 = fftshift(gyro_y_50);
% n_50 = length(gyro_y_50);                    % number of samples
% f_50 = (-n_50/2:(n_50/2)-1)*(DS.fs/n_50);    % frequency range
% y_power_50 = abs(gyro_y_50).^2/n_50; 
% 
% %60
% gyro_y_60 = fft(DS.gyro_y(i_ini_60:i_fin_60)-mean(DS.gyro_y(i_ini_60:i_fin_60)));
% gyro_y_60 = fftshift(gyro_y_60);
% n_60 = length(gyro_y_60);                    % number of samples
% f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
% y_power_60 = abs(gyro_y_60).^2/n_60; 
% 
% 
% 
% %70
% gyro_y_70 = fft(DS.gyro_y(i_ini_70:i_fin_70)-mean(DS.gyro_y(i_ini_70:i_fin_70)));
% gyro_y_70 = fftshift(gyro_y_70);
% n_70 = length(gyro_y_70);                    % number of samples
% f_70 = (-n_70/2:(n_70/2)-1)*(DS.fs/n_70);    % frequency range
% y_power_70 = abs(gyro_y_70).^2/n_70;  
% 
% %80
% gyro_y_80 = fft(DS.gyro_y(i_ini_80:i_fin_80)-mean(DS.gyro_y(i_ini_80:i_fin_80)));
% gyro_y_80 = fftshift(gyro_y_80);
% n_80 = length(gyro_y_80);                    % number of samples
% f_80 = (-n_80/2:(n_80/2)-1)*(DS.fs/n_80);    % frequency range
% y_power_80 = abs(gyro_y_80).^2/n_80; 
% 
% %90
% gyro_y_90 = fft(DS.gyro_y(i_ini_90:i_fin_90)-mean(DS.gyro_y(i_ini_90:i_fin_90)));
% gyro_y_90 = fftshift(gyro_y_90);
% n_90 = length(gyro_y_90);                    % number of samples
% f_90 = (-n_90/2:(n_90/2)-1)*(DS.fs/n_90);    % frequency range
% y_power_90 = abs(gyro_y_90).^2/n_90; 
% 
% 
% % %% Powers gyro z
% % 
% %50
% 
% gyro_z_50 = fft(DS.gyro_z(i_ini_50:i_fin_50)-mean(DS.gyro_z(i_ini_50:i_fin_50)));
% gyro_z_50 = fftshift(gyro_z_50);
% n_50 = length(gyro_z_50);                    % number of samples
% f_50 = (-n_50/2:(n_50/2)-1)*(DS.fs/n_50);    % frequency range
% z_power_50 = abs(gyro_z_50).^2/n_50; 
% 
% %60
% gyro_z_60 = fft(DS.gyro_z(i_ini_60:i_fin_60)-mean(DS.gyro_z(i_ini_60:i_fin_60)));
% gyro_z_60 = fftshift(gyro_z_60);
% n_60 = length(gyro_z_60);                    % number of samples
% f_60 = (-n_60/2:(n_60/2)-1)*(DS.fs/n_60);    % frequency range
% z_power_60 = abs(gyro_z_60).^2/n_60; 
% 
% 
% 
% %70
% gyro_z_70 = fft(DS.gyro_z(i_ini_70:i_fin_70)-mean(DS.gyro_z(i_ini_70:i_fin_70)));
% gyro_z_70 = fftshift(gyro_z_70);
% n_70 = length(gyro_z_70);                    % number of samples
% f_70 = (-n_70/2:(n_70/2)-1)*(DS.fs/n_70);    % frequency range
% z_power_70 = abs(gyro_z_70).^2/n_70;  
% 
% %80
% gyro_z_80 = fft(DS.gyro_z(i_ini_80:i_fin_80)-mean(DS.gyro_z(i_ini_80:i_fin_80)));
% gyro_z_80 = fftshift(gyro_z_80);
% n_80 = length(gyro_z_80);                    % number of samples
% f_80 = (-n_80/2:(n_80/2)-1)*(DS.fs/n_80);    % frequency range
% z_power_80 = abs(gyro_z_80).^2/n_80; 
% 
% %90
% gyro_z_90 = fft(DS.gyro_z(i_ini_90:i_fin_90)-mean(DS.gyro_z(i_ini_90:i_fin_90)));
% gyro_z_90 = fftshift(gyro_z_90);
% n_90 = length(gyro_z_90);                    % number of samples
% f_90 = (-n_90/2:(n_90/2)-1)*(DS.fs/n_90);    % frequency range
% z_power_90 = abs(gyro_z_90).^2/n_90; 
% 
% 
% 
% 
% % Total freq power 
% tot_power_50 = x_power_50+y_power_50+z_power_50;
% tot_power_60 = x_power_60+y_power_60+z_power_60;
% tot_power_70 = x_power_70+y_power_70+z_power_70;
% tot_power_80 = x_power_80+y_power_80+z_power_80;
% tot_power_90 = x_power_90+y_power_90+z_power_90;
% 
% % 
% % 
% % Plotting power;
% 
% figure
% subplot(5,1,1)
% plot(f_50,tot_power_50,'-b');grid on;
% xlabel('freq'); ylabel('power');
% axis([0 5 0 0.2])
% subplot(5,1,2)
% plot(f_60,tot_power_60,'-r');grid on;
% xlabel('freq'); ylabel('power');
% axis([0 5 0 0.2])
% subplot(5,1,3)
% plot(f_70,tot_power_70,'-k');grid on;
% xlabel('freq'); ylabel('power');
% axis([0 5 0 0.2])
% subplot(5,1,4)
% plot(f_80,tot_power_80,'-m');grid on;
% xlabel('freq'); ylabel('power');
% axis([0 5 0 0.2])
% subplot(5,1,5)
% plot(f_90,tot_power_90,'-g');grid on;
% xlabel('freq'); ylabel('power');
% axis([0 5 0 0.2])
% % saveas(gca,'frequency_power_per','svg')
% %Plotting the gyroscope
% 
% % 
% % figure
% % subplot(3,1,1);
% % plot(DS.time(i_ini_50:i_fin_50),DS.gyro_x(i_ini_50:i_fin_50),'-b');
% % grid on; hold on;
% % plot(DS.time(i_ini_60:i_fin_60),DS.gyro_x(i_ini_60:i_fin_60),'-r');
% % plot(DS.time(i_ini_70:i_fin_70),DS.gyro_x(i_ini_70:i_fin_70),'-k');
% % plot(DS.time(i_ini_80:i_fin_80),DS.gyro_x(i_ini_80:i_fin_80),'-m');
% % plot(DS.time(i_ini_90:i_fin_90),DS.gyro_x(i_ini_90:i_fin_90),'-g');
% % xlabel('temps');
% % ylabel('gyrox');
% % hold on;
% % axis([DS.time(i_ini_50) DS.time(i_fin_90) min(DS.gyro_x) max(DS.gyro_x)]);
% % % for i =1:18
% % %     plot([10*(i-1)+t_ini,10*(i-1)+t_ini],[-5,5],'-r');
% % % end
% % % axis([0 220 -5 5]);
% % 
% % 
% % subplot(3,1,2)
% % plot(DS.time(i_ini_50:i_fin_50),DS.gyro_y(i_ini_50:i_fin_50),'-b');
% % grid on; hold on;
% % plot(DS.time(i_ini_60:i_fin_60),DS.gyro_y(i_ini_60:i_fin_60),'-r');
% % plot(DS.time(i_ini_70:i_fin_70),DS.gyro_y(i_ini_70:i_fin_70),'-k');
% % plot(DS.time(i_ini_80:i_fin_80),DS.gyro_y(i_ini_80:i_fin_80),'-m');
% % plot(DS.time(i_ini_90:i_fin_90),DS.gyro_y(i_ini_90:i_fin_90),'-g');
% % xlabel('temps');
% % ylabel('gyroy');
% % hold on;
% % axis([DS.time(i_ini_50) DS.time(i_fin_90) min(DS.gyro_y) max(DS.gyro_y)]);
% % % for i =1:18
% % %     plot([10*(i-1)+t_ini,10*(i-1)+t_ini],[-5,5],'-r');
% % % end
% % % axis([DS.time(i_ini_50) DS.time(i_fin_90) ...
% % %      min(DS.gyro_y(i_ini_50:i_fin_50) 5]);
% % 
% % 
% % subplot(3,1,3)
% % plot(DS.time(i_ini_50:i_fin_50),DS.gyro_z(i_ini_50:i_fin_50),'-b');
% % grid on; hold on;
% % plot(DS.time(i_ini_60:i_fin_60),DS.gyro_z(i_ini_60:i_fin_60),'-r');
% % plot(DS.time(i_ini_70:i_fin_70),DS.gyro_z(i_ini_70:i_fin_70),'-k');
% % plot(DS.time(i_ini_80:i_fin_80),DS.gyro_z(i_ini_80:i_fin_80),'-m');
% % plot(DS.time(i_ini_90:i_fin_90),DS.gyro_z(i_ini_90:i_fin_90),'-g');
% % xlabel('temps');
% % ylabel('gyroz');
% % hold on;
% % axis([DS.time(i_ini_50) DS.time(i_fin_90) min(DS.gyro_z) max(DS.gyro_z)]);
% % % for i =1:18
% % %     plot([10*(i-1)+t_ini,10*(i-1)+t_ini],[-5,5],'-r');
% % % end
% % % axis([0 220 -5 5]);
% 
% % saveas(gca,'gyro_xyz','svg');
% 
% 
% 
% 
% % % Plotting the acceleration
% % figure
% % subplot(3,1,1)
% % plot(DS.time,DS.accel_xf,'-b'); grid on;
% % xlabel('temps');
% % ylabel('accelx');
% % hold on;
% % for i =1:18
% %     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% % end
% % 
% % 
% % subplot(3,1,2)
% % plot(DS.time,DS.accel_yf,'-b'); grid on;
% % xlabel('temps');
% % ylabel('accely');
% % hold on;
% % for i =1:18
% %     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% % end
% % 
% % 
% % subplot(3,1,3)
% % plot(DS.time,DS.accel_zf,'-b'); grid on;
% % xlabel('temps');
% % ylabel('accelz');
% % hold on;
% % for i =1:18
% %     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% % end
% % 
% % 
% % 
% % % Plotting the orientation
% % figure
% % subplot(3,1,1)
% % plot(DS.time,DS.roll,'-b');grid on;
% % xlabel('temps');
% % ylabel('roll');
% % hold on;
% % for i =1:18
% %     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% % end
% % subplot(3,1,2)
% % plot(DS.time,DS.pitch,'-b');grid on;
% % xlabel('temps');
% % ylabel('pitch');
% % hold on;
% % for i =1:18
% %     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% % end
% % subplot(3,1,3)
% % plot(DS.time,DS.yaw,'-b');grid on;
% % xlabel('temps');
% % ylabel('yaw');
% % hold on;
% % for i =1:18
% %     plot([10*(i-1)+35,10*(i-1)+35],[-5,5],'-r');
% % end
% % 
% % 
% % 
% % 
