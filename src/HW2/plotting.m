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
% h = csvread('../../data/
OEdrag = csvread('../../data/OEhistDrag_HW2.csv');
Edrag = csvread('../../data/EhistDrag_HW2.csv');

% change in OE
dOE = OE - OEdrag;

%plot OEs
figure
subplot(3,2,1)
plot(t/3600,OE(1,:),'LineWidth',2)
ylabel('$a$ [km]')
% xlabel('Time [hours]')
grid on
% saveas(gcf,[savepath 'semimajor.pdf'])

subplot(3,2,2)
plot(t/3600, OE(2,:), 'LineWidth', 2)
% xlabel('Time [hours]')
ylabel('$e$')
grid on
% saveas(gcf,[savepath 'eccentricity.pdf'])

subplot(3,2,3)
plot(t/3600, OE(3,:), 'LineWidth', 2)
% xlabel('Time [hours]')
ylabel('$i$ [deg]')
grid on
% saveas(gcf,[savepath 'inclination.pdf'])

subplot(3,2,4)
plot(t/3600, OE(4,:), 'LineWidth', 2)
% xlabel('Time [hours]')
ylabel('$\Omega$ [deg]')
grid on
% saveas(gcf,[savepath 'LongAscend .pdf'])

subplot(3,2,5)
plot(t/3600, OE(5,:), 'LineWidth', 2)
xlabel('Time [hours]')
grid on
ylabel('$\omega$ [deg]')

subplot(3,2,6)
plot(t/3600, OE(6,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$T_p$ [sec]')
grid on
saveas(gcf,[savepath 'OE.pdf'])

figure
plot(t/3600, OE(7,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$P$ [sec]')
grid on
saveas(gcf,[savepath 'OrbitPeriod.pdf'])

%plot energy
figure
plot(t/3600, E(1,:) - E(2,:) - E(1,1) + E(2,1), 'LineWidth',2)
grid on
xlabel('Time [hours]')
ylabel('Change in Specific Energy [J/kg]')
saveas(gcf,[savepath 'energy.pdf'])

%plot change in OEs
figure
subplot(3,2,1)
plot(t/3600,dOE(1,:),'LineWidth',2)
ylabel('$a$ [km]')
% xlabel('Time [hours]')
grid on
saveas(gcf,[savepath 'semimajor.pdf'])

subplot(3,2,2)
plot(t/3600, dOE(2,:), 'LineWidth', 2)
% xlabel('Time [hours]')
ylabel('$e$')
grid on
% saveas(gcf,[savepath 'eccentricity.pdf'])

subplot(3,2,3)
plot(t/3600, dOE(3,:), 'LineWidth', 2)
% xlabel('Time [hours]')
ylabel('$i$ [deg]')
grid on
% saveas(gcf,[savepath 'inclination.pdf'])

subplot(3,2,4)
plot(t/3600, dOE(4,:), 'LineWidth', 2)
% xlabel('Time [hours]')
ylabel('$\Omega$ [deg]')
grid on
% saveas(gcf,[savepath 'LongAscend .pdf'])

subplot(3,2,5)
plot(t/3600, dOE(5,:), 'LineWidth', 2)
xlabel('Time [hours]')
grid on
ylabel('$\omega$ [deg]')

subplot(3,2,6)
plot(t/3600, dOE(6,:), 'LineWidth', 2)
xlabel('Time [hours]')
ylabel('$T_p$ [sec]')
grid on
saveas(gcf,[savepath 'dOE.pdf'])

%plot energy with drag
figure
plot(t/3600, Edrag(1,:) - Edrag(2,:) - Edrag(1,1) + Edrag(2,1), 'LineWidth',2)
grid on
xlabel('Time [hours]')
ylabel('Change in Specific Energy [J/kg]')
saveas(gcf,[savepath 'energydrag.pdf'])

