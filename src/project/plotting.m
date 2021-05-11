clear
close all
clc


%latex
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%savepath for images
savepath = '../../doc/project/figs/';

%matlab scripts
addpath('../../../matlabScripts/');

%load data
xhat = csvread("../../data/xhat_proj.csv");
Phat = csvread("../../data/Phat_proj.csv");
prefit = csvread("../../data/prefit_res_proj.csv");
postfit = csvread("../../data/postfit_res_proj.csv");
Pyy_pre = csvread("../../data/Pyy_pre_proj.csv");
Pyy_post = csvread("../../data/Pyy_post_proj.csv");
z = csvread("../../data/meas_proj_set1.csv");
% z = csvread("../../data/meas_proj_set3.csv");
xhatA = csvread("../../data/xhat_A_NAG.csv");
PhatA = csvread("../../data/Phat_A_NAG.csv");
xhatB = csvread("../../data/xhat_B_NAG.csv");
PhatB = csvread("../../data/Phat_B_NAG.csv");
xhatC = csvread("../../data/xhat_C_NAG.csv");
PhatC = csvread("../../data/Phat_C_NAG.csv");
xhatD = csvread("../../data/xhat_D_NAG.csv");
PhatD = csvread("../../data/Phat_D_NAG.csv");
xhatE = csvread("../../data/xhat_E_NAG.csv");
PhatE = csvread("../../data/Phat_E_NAG.csv");
xhatF = csvread("../../data/xhat_F_NAG.csv");
PhatF = csvread("../../data/Phat_F_NAG.csv");
xhatG = csvread("../../data/xhat_G_NAG.csv");
PhatG = csvread("../../data/Phat_G_NAG.csv");

%data for NAG
marcus_pos_caseA = xhatA(1:3);
marcus_pos_caseB = xhatB(1:3);
marcus_pos_caseC = xhatC(1:3);
marcus_pos_caseD = xhatD(1:3);
marcus_pos_caseE = xhatE(1:3);
marcus_pos_caseF = xhatF(1:3);
marcus_pos_caseG = xhatG(1:3);
marcus_poscov_caseA = PhatA(1:3,1:3);
marcus_poscov_caseB = PhatB(1:3,1:3);
marcus_poscov_caseC = PhatC(1:3,1:3);
marcus_poscov_caseD = PhatD(1:3,1:3);
marcus_poscov_caseE = PhatE(1:3,1:3);
marcus_poscov_caseF = PhatF(1:3,1:3);
marcus_poscov_caseG = PhatG(1:3,1:3);
save("marcus.mat","marcus*");

%target indexes for plotting
% idxs = 100:length(z(:,2));
idxs = 1:length(z(:,2));

%true position velocity after 24 hours
truth24hours = [-6330.16736001325;
    3306.81591178162;
    127.736863438565;
    -3.43787907681733;
    -6.63350511630163;
    -0.235613730204275];

%get orbital elements for best estimate
best_est = truth24hours;
% best_est = xhatF;
d2r = pi/180;
[~, ~, i, Ohm, w, theta] = FunState2OE(best_est(1:3),best_est(4:6));

%get rotation between ECI and RSW
R_ECI2PQW = angle2dcm(Ohm*d2r,i*d2r,w*d2r,'ZXZ');
R_PQW2RSW = angle2dcm(-theta*d2r,0,0,'ZYX');
R_total = R_PQW2RSW*R_ECI2PQW;

%find deviations from best estimate
devA = R_total*(xhatA(1:3) - best_est(1:3));
devB = R_total*(xhatB(1:3) - best_est(1:3));
devC = R_total*(xhatC(1:3) - best_est(1:3));
devD = R_total*(xhatD(1:3) - best_est(1:3));
devE = R_total*(xhatE(1:3) - best_est(1:3));
devF = R_total*(xhatF(1:3) - best_est(1:3));
% devG = R_total*(xhatG(1:3) - best_est(1:3));
devG = 0*R_total*(xhatG(1:3) - best_est(1:3));
PhatdevA = R_total*PhatA*R_total';
PhatdevB = R_total*PhatB*R_total';
PhatdevC = R_total*PhatC*R_total';
PhatdevD = R_total*PhatD*R_total';
PhatdevE = R_total*PhatE*R_total';
PhatdevF = R_total*PhatF*R_total';
% PhatdevG = R_total*PhatG*R_total';
PhatdevG = 0*R_total*PhatG*R_total';

% [~, ~, i, Ohm, w, theta] = FunState2OE([100; .001; .001],[0; 100; 0]);
% R_ECI2PQW = angle2dcm(Ohm*d2r,i*d2r,w*d2r,'ZXZ');
% R_PQW2RSW = angle2dcm(-theta*d2r,0,0,'ZYX');
% R_total = R_PQW2RSW*R_ECI2PQW;
% R_total*[0; 0; 1]

%measurement standard deviations
sqrtR = zeros(2,length(z(:,1)));
station1idx = find(z(:,1) == 1);
station2idx = find(z(:,1) == 2);
station3idx = find(z(:,1) == 3);
sqrtR(1,station1idx) = 10/1000;
sqrtR(1,station2idx) = 5/1000;
sqrtR(1,station3idx) = 10/1000;
sqrtR(2,station1idx) = 0.5/(1000*1000);
sqrtR(2,station2idx) = 1/(1000*1000);
sqrtR(2,station3idx) = 0.5/(1000*1000);

%plot residuals

figure
subplot(2,1,1)
scatter(z(station1idx,2)/3600,1000*postfit(1,station1idx),4,'filled')
hold on
scatter(z(station2idx,2)/3600,1000*postfit(1,station2idx),4,'filled')
scatter(z(station3idx,2)/3600,1000*postfit(1,station3idx),4,'filled')
plot(z(idxs,2)/3600,1000*3*sqrt(Pyy_post(1,idxs)),'r')
plot(z(idxs,2)/3600,1000*-3*sqrt(Pyy_post(1,idxs)),'r')
title('Post-Fit Residuals')
grid on
ylabel('Range [m]')
legend('Station 1','Station 2','Station 3')

subplot(2,1,2)
scatter(z(station1idx,2)/3600,100*1000*postfit(2,station1idx),4,'filled')
hold on
scatter(z(station2idx,2)/3600,100*1000*postfit(2,station2idx),4,'filled')
scatter(z(station3idx,2)/3600,100*1000*postfit(2,station3idx),4,'filled')
plot(z(idxs,2)/3600,100*1000*3*sqrt(Pyy_post(4,idxs)),'r')
plot(z(idxs,2)/3600,100*1000*-3*sqrt(Pyy_post(4,idxs)),'r')
ylabel('Range Rate [mm/sec]')
xlabel('Time (hours)')
grid on
legend('Station 1','Station 2','Station 3')

figure
subplot(2,1,1)
scatter(z(station1idx,2)/3600,1000*prefit(1,station1idx),4,'filled')
hold on
scatter(z(station2idx,2)/3600,1000*prefit(1,station2idx),4,'filled')
scatter(z(station3idx,2)/3600,1000*prefit(1,station3idx),4,'filled')
plot(z(idxs,2)/3600,1000*3*sqrt(Pyy_pre(1,idxs)),'r')
plot(z(idxs,2)/3600,1000*-3*sqrt(Pyy_pre(1,idxs)),'r')
title('Pre-Fit Residuals')
grid on
ylabel('Range [m]')
legend('Station 1','Station 2','Station 3')

subplot(2,1,2)
scatter(z(station1idx,2)/3600,100*1000*prefit(2,station1idx),4,'filled')
hold on
scatter(z(station2idx,2)/3600,100*1000*prefit(2,station2idx),4,'filled')
scatter(z(station3idx,2)/3600,100*1000*prefit(2,station3idx),4,'filled')
plot(z(idxs,2)/3600,100*1000*3*sqrt(Pyy_pre(4,idxs)),'r')
plot(z(idxs,2)/3600,100*1000*-3*sqrt(Pyy_pre(4,idxs)),'r')
grid on
ylabel('Range Rate [mm/sec]')
xlabel('Time (hours)')
legend('Station 1','Station 2','Station 3')

figure
subplot(2,1,1)
scatter(z(station1idx,2)/3600,1000*(postfit(1,station1idx) - prefit(1,station1idx)),4,'filled')
hold on
scatter(z(station2idx,2)/3600,1000*(postfit(1,station2idx) - prefit(1,station2idx)),4,'filled')
scatter(z(station3idx,2)/3600,1000*(postfit(1,station3idx) - prefit(1,station3idx)),4,'filled')
title('Update in Residuals')
grid on
ylabel('Range [m]')
legend('Station 1','Station 2','Station 3')

subplot(2,1,2)
scatter(z(station1idx,2)/3600,100*1000*(postfit(2,station1idx) - prefit(2,station1idx)),4,'filled')
hold on
scatter(z(station2idx,2)/3600,100*1000*(postfit(2,station2idx) - prefit(2,station2idx)),4,'filled')
scatter(z(station3idx,2)/3600,100*1000*(postfit(2,station3idx) - prefit(2,station3idx)),4,'filled')
grid on
ylabel('Range Rate [mm/sec]')
xlabel('Time (hours)')
legend('Station 1','Station 2','Station 3')

% figure
% plot3(xhat(1,:),xhat(2,:),xhat(3,:))
% title('Position Estimate')
% xlabel('x')
% ylabel('y')
% zlabel('z')

% figure
% plot(z(:,2),xhat(3,:))
% title('Position Estimate')

%plot estimates
xyfig = figure;
sigval = 3;
scatter(devA(1),devA(2),'x')
hold on
scatter(devB(1),devB(2),'x')
scatter(devC(1),devC(2),'x')
scatter(devD(1),devD(2),'x')
scatter(devE(1),devE(2),'x')
scatter(devF(1),devF(2),'x')
scatter(devG(1),devG(2),'x')
plot_elipse(xyfig,PhatdevA([1 2],[1 2]),devA([1 2]),sigval,'',false)
plot_elipse(xyfig,PhatdevB([1 2],[1 2]),devB([1 2]),sigval,'',false)
plot_elipse(xyfig,PhatdevC([1 2],[1 2]),devC([1 2]),sigval,'',false)
plot_elipse(xyfig,PhatdevD([1 2],[1 2]),devD([1 2]),sigval,'',false)
plot_elipse(xyfig,PhatdevE([1 2],[1 2]),devE([1 2]),sigval,'',false)
plot_elipse(xyfig,PhatdevF([1 2],[1 2]),devF([1 2]),sigval,'',false)
plot_elipse(xyfig,PhatdevG([1 2],[1 2]),devG([1 2]),sigval,'',false)
xlabel('Radial [km]')
ylabel('In Track [km]')
text(devA(1)+.005,devA(2),'A')
text(devB(1)+.005,devB(2)-.005,'B')
text(devC(1)+.005,devC(2),'C')
text(devD(1)+.005,devD(2),'D')
text(devE(1)-.005,devE(2),'E')
text(devF(1)+.005,devF(2),'F')
text(devG(1)+.005,devG(2)+.005,'G')

xzfig = figure;
scatter(devA(1),devA(3),'x')
hold on
scatter(devB(1),devB(3),'x')
scatter(devC(1),devC(3),'x')
scatter(devD(1),devD(3),'x')
scatter(devE(1),devE(3),'x')
scatter(devF(1),devF(3),'x')
scatter(devG(1),devG(3),'x')
plot_elipse(xzfig,PhatdevA([1 3],[1 3]),devA([1 3]),sigval,'',false)
plot_elipse(xzfig,PhatdevB([1 3],[1 3]),devB([1 3]),sigval,'',false)
plot_elipse(xzfig,PhatdevC([1 3],[1 3]),devC([1 3]),sigval,'',false)
plot_elipse(xzfig,PhatdevD([1 3],[1 3]),devD([1 3]),sigval,'',false)
plot_elipse(xzfig,PhatdevE([1 3],[1 3]),devE([1 3]),sigval,'',false)
plot_elipse(xzfig,PhatdevF([1 3],[1 3]),devF([1 3]),sigval,'',false)
plot_elipse(xzfig,PhatdevG([1 3],[1 3]),devG([1 3]),sigval,'',false)
xlabel('Radial [km]')
ylabel('Cross Track [km]')
text(devA(1)+.005,devA(3),'A')
text(devB(1)+.005,devB(3)-.005,'B')
text(devC(1)+.005,devC(3),'C')
text(devD(1)+.005,devD(3),'D')
text(devE(1)-.005,devE(3),'E')
text(devF(1)+.005,devF(3),'F')
text(devG(1)+.005,devG(3)+.005,'G')

yzfig = figure;
scatter(devA(2),devA(3),'x')
hold on
scatter(devB(2),devB(3),'x')
scatter(devC(2),devC(3),'x')
scatter(devD(2),devD(3),'x')
scatter(devE(2),devE(3),'x')
scatter(devF(2),devF(3),'x')
scatter(devG(2),devG(3),'x')
plot_elipse(yzfig,PhatdevA([2 3],[2 3]),devA([2 3]),sigval,'',false)
plot_elipse(yzfig,PhatdevB([2 3],[2 3]),devB([2 3]),sigval,'',false)
plot_elipse(yzfig,PhatdevC([2 3],[2 3]),devC([2 3]),sigval,'',false)
plot_elipse(yzfig,PhatdevD([2 3],[2 3]),devD([2 3]),sigval,'',false)
plot_elipse(yzfig,PhatdevE([2 3],[2 3]),devE([2 3]),sigval,'',false)
plot_elipse(yzfig,PhatdevF([2 3],[2 3]),devF([2 3]),sigval,'',false)
plot_elipse(yzfig,PhatdevG([2 3],[2 3]),devG([2 3]),sigval,'',false)
xlabel('In Track [km]')
ylabel('Cross Track [km]')
text(devA(2)+.005,devA(3),'A')
text(devB(2)+.005,devB(3)-.005,'B')
text(devC(2)+.005,devC(3),'C')
text(devD(2)+.005,devD(3),'D')
text(devE(2)-.005,devE(3),'E')
text(devF(2)+.005,devF(3),'F')
text(devG(2)+.005,devG(3)+.005,'G')

