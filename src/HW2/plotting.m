%This script does all the plotting of my data, because plotting in cpp is hard :/

%startup
clear
close all
clc

%latex
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%savepath for images
savepath = '../../doc/HW2/figs/';

%load the data
OE = csvread('../../data/OEhist_HW2.csv');
t = csvread('../../data/thist_HW2.csv');
E = csvread('../../data/Ehist_HW2.csv');

%plot OEs
figure
plot(t/3600,OE(1,:),'LineWidth',2)
ylabel('$a$ [km]')
xlabel('Time [hours]')
grid on
saveas(gcf,[savepath 'semimajor.pdf'])

figure
plot(t/3600, OE(2,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$e$')
saveas(gcf,[savepath 'eccentricity.pdf'])

figure
plot(t/3600, OE(3,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$i$ [deg]')
saveas(gcf,[savepath 'inclination.pdf'])

figure
plot(t/3600, OE(4,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$\Omega$ [deg]')
saveas(gcf,[savepath 'LongAscend .pdf'])

figure
plot(t/3600, OE(5,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$\omega$ [deg]')
saveas(gcf,[savepath 'ArgPeriaps.pdf'])

figure
plot(t/3600, OE(6,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$T_p$ [sec]')
saveas(gcf,[savepath 'semimajor.pdf'])

figure
plot(t/3600, OE(7,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$P$ [sec]')
