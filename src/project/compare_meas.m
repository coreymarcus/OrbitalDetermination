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

%plot
figure
subplot(2,1,1)
scatter(t, LEO_DATA_Apparent(:,3));
hold on
scatter(t,pred_meas(1,:))


figure
subplot(2,1,1)
scatter(t,err(1,:))

subplot(2,1,2)
scatter(t,err(2,:))