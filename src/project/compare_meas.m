clear
close all
clc

%load data
pred_meas = csvread("../../data/pred_meas.csv");
load('../../data/LEO_DATA_Apparent.mat');

%error
err = LEO_DATA_Apparent(:,3:4)' - pred_meas;

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