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
Pyy = csvread("../../data/Pyy_proj.csv");
z = csvread("../../data/meas_proj_set1.csv");

%target indexes for plotting
idxs = 100:length(z(:,2));

%plot residuals

figure
subplot(2,1,1)
scatter(z(idxs,2),postfit(1,idxs),4,'filled')
hold on
plot(z(idxs,2),3*sqrt(Pyy(1,idxs)),'r')
plot(z(idxs,2),-3*sqrt(Pyy(1,idxs)),'r')
title('Post-Fit Residuals')
grid on

subplot(2,1,2)
scatter(z(idxs,2),postfit(2,idxs),4,'filled')
grid on
hold on
plot(z(idxs,2),3*sqrt(Pyy(4,idxs)),'r')
plot(z(idxs,2),-3*sqrt(Pyy(4,idxs)),'r')

figure
subplot(2,1,1)
scatter(z(idxs,2),prefit(1,idxs),4,'filled')
title('Pre-Fit Residuals')
grid on

subplot(2,1,2)
scatter(z(idxs,2),prefit(2,idxs),4,'filled')
grid on

% figure
% plot3(xhat(1,:),xhat(2,:),xhat(3,:))
% title('Position Estimate')
% xlabel('x')
% ylabel('y')
% zlabel('z')

% figure
% plot(z(:,2),xhat(3,:))
% title('Position Estimate')
