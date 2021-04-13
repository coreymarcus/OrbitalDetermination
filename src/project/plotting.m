clear
close all
clc


%latex
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%savepath for images
savepath = '../../doc/project/figs/';

%load data
xhat = csvread("../../data/xhat_proj.csv");
Phat = csvread("../../data/Phat_proj.csv");
prefit = csvread("../../data/prefit_res_proj.csv");
postfit = csvread("../../data/postfit_res_proj.csv");
z = csvread("../../data/meas_proj_set1.csv");

%plot residuals

figure
subplot(2,1,1)
scatter(z(:,2),postfit(1,:),4,'filled')
title('Post-Fit Residuals')

subplot(2,1,2)
scatter(z(:,2),postfit(2,:),4,'filled')

figure
subplot(2,1,1)
scatter(z(:,2),prefit(1,:),4,'filled')
title('Pre-Fit Residuals')

subplot(2,1,2)
scatter(z(:,2),prefit(2,:),4,'filled')

% figure
% plot3(xhat(1,:),xhat(2,:),xhat(3,:))
% title('Position Estimate')
% xlabel('x')
% ylabel('y')
% zlabel('z')

figure
plot(z(:,2),xhat(3,:))
title('Position Estimate')
