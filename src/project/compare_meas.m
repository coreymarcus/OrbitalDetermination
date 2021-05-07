clear
close all
clc

%load data
pred_meas = csvread("../../data/pred_meas.csv");
load('../../data/LEO_DATA_Apparent.mat');

%load data for ECI stuff
load('../../data/R_test.mat','r_ECI');
r_ECI1 = r_ECI;
clear r_ECI;
load('../../data/connerRotations.mat');
r_ECI2 = rECI;
clear rECI;
pred_ECI = csvread("../../data/pred_ECIpos.csv")';

%error
err = LEO_DATA_Apparent(:,3:4)' - pred_meas;

%ECI error
err1 = 100*1000*(r_ECI1 - pred_ECI);
err2 = 100*1000*(r_ECI2 - pred_ECI);
err3 = 100*1000*(r_ECI1 - r_ECI2);

%time
t = LEO_DATA_Apparent(:,2);

%station idx
station1idx = find(LEO_DATA_Apparent(:,1) == 1);
station2idx = find(LEO_DATA_Apparent(:,1) == 2);
station3idx = find(LEO_DATA_Apparent(:,1) == 3);

%plot
% figure
% subplot(2,1,1)
% hold on
% scatter(t(station1idx), LEO_DATA_Apparent(station1idx,3));
% scatter(t(station2idx), LEO_DATA_Apparent(station2idx,3));
% scatter(t(station3idx), LEO_DATA_Apparent(station3idx,3));
% scatter(t,pred_meas(1,:))


figure
subplot(2,1,1)
hold on
scatter(t(station1idx),1000*err(1,station1idx),8,'filled')
scatter(t(station2idx),1000*err(1,station2idx),8,'filled')
scatter(t(station3idx),1000*err(1,station3idx),8,'filled')
title('Pred Measurement Error')
ylabel('Range [m]')

subplot(2,1,2)
hold on
scatter(t(station1idx),100*1000*err(2,station1idx),8,'filled')
scatter(t(station2idx),100*1000*err(2,station2idx),8,'filled')
scatter(t(station3idx),100*1000*err(2,station3idx),8,'filled')
ylabel('Range-Rate [mm/sec]')
xlabel('Time [sec]')

figure
subplot(3,1,1)
hold on
plot(err1(:,1))
plot(err2(:,1))
ylabel('Error (cm)');
title("ECEF2ECI Comparison with Other Student")

subplot(3,1,2)
hold on
plot(err1(:,2))
plot(err2(:,2))
ylabel('Error (cm)');

subplot(3,1,3)
hold on
plot(err1(:,3))
plot(err2(:,3))
ylabel('Error (cm)');

% figure
% subplot(3,1,1)
% hold on
% plot(err3(:,1))
% 
% subplot(3,1,2)
% hold on
% plot(err3(:,2))
% 
% subplot(3,1,3)
% hold on
% plot(err3(:,3))