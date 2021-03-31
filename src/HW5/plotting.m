clear
close all
clc

%latex
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%savepath for images
savepath = '../../doc/HW5/figs/';

%load data
A0 = csvread("../../data/myA0_HW5.csv");
H0 = csvread("../../data/myH0_HW5.csv");
Phi = csvread("../../data/myPhi_HW5.csv");
zbar = csvread("../../data/mymeas_HW5.csv");

A0_true = load('A_t0.mat');
H0_true = load('H_Tilde_t0.mat');
Phi_true = load('Phi_21600_0.mat');
z_true = load('LEO_DATA_Apparent.mat');


%remove structure
A0_true = A0_true.A;
H0_true = H0_true.H_TILDA;
Phi_true = Phi_true.PHI_t_120;
z_true = z_true.LEO_DATA_Apparent;

%Problem 1
Hreldiff = abs((H0 - H0_true(:,1:6))./H0_true(:,1:6));
Areldiff = abs((A0 - A0_true(1:6,1:6))./A0_true(1:6,1:6));

Hreldiff(isnan(Hreldiff)) = 0;
Areldiff(isnan(Areldiff)) = 0;

figure
histogram(reshape(log10(Areldiff),36,1),5)
title('Histogram of log10 of relative difference in A')
xlabel('$log10((A - A_{true})./A_{true})$')
ylabel('Frequency')

%problem 2
Phireldiff = abs((Phi - Phi_true(1:6,1:6))./Phi_true(1:6,1:6));
Phireldiff(isnan(Phireldiff)) = 0;


figure
histogram(reshape(log10(Phireldiff),36,1),5)
title('Histogram of log10 of relative difference in $\Phi$')
xlabel('$log10((\Phi - \Phi_{true})./\Phi_{true})$')
ylabel('Frequency')

% problem 3

%residual
res = zbar' - z_true(:,3:4);

%time
t = z_true(:,2)/3600;

figure
subplot(2,1,1)
scatter(t,res(:,1),'*')
title('Pre-Fit Measurement Residuals')
ylabel('Range Residual (km)')

subplot(2,1,2)
scatter(t,res(:,2),'*');
xlabel('Time (hours)')
ylabel('Range Rate Residual (km/sec)')
