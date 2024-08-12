clear all
close all
clc
load Meas_Dist_all.mat
Temp_Dist=Meas_Dist_all(:,1);
Rad_Dist=Meas_Dist_all(:,2);
total_sim_time=20*60*60-1; %20 hours
ts=60;
Tm=ts;
tvec=(0:ts:total_sim_time)';
Temp_vec=[tvec Temp_Dist(1:(total_sim_time/ts)+1)];
Rad_vec=[tvec Rad_Dist(1:(total_sim_time/ts)+1)];
pH_sp=8*ones(length(tvec),1); % Define Setpoint

% pH_sp(1:400)=7.2;
% pH_sp(850:end)=7.5;

plantdef_reactor;
open('reactor_cl_test_21a.slx');
sim('reactor_cl_test_21a.slx');




figure;
subplot(4,1,1);plot(tvec/(60*60),pH_out);ylabel('pH');hold on;plot(tvec/(60*60),pH_sp);subplot(4,1,2);plot(tvec/(60*60),CO2_in);ylabel('CO_2 (L/min)');subplot(4,1,3);plot(Temp_vec(:,1)/(60*60),Temp_vec(:,2));ylabel('Temp (C)');subplot(4,1,4);plot(Rad_vec(:,1)/(60*60),Rad_vec(:,2));ylabel('Rad (W/m^2)');xlabel('Time (hours)')