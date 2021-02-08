%This script does all the plotting of my data, because plotting in cpp is hard :/

%startup
clear
close all
clc

%load the x history
xhist = csvread('../../data/xhist_HW1.csv');

%plot position and velocity and accel
figure
plot3(xhist(1,:), xhist(2,:), xhist(3,:))
title('Position History','Interpreter','latex')
xlabel('x [km]','Interpreter','latex')
ylabel('y [km]','Interpreter','latex')
zlabel('z [km]','Interpreter','latex')
grid on

figure
plot3(xhist(4,:), xhist(5,:), xhist(6,:))
title('Velocity History','Interpreter','latex')
xlabel('x [km/sec]','Interpreter','latex')
ylabel('y [km/sec]','Interpreter','latex')
zlabel('z [km/sec]','Interpreter','latex')
grid on

figure
plot3(xhist(7,:), xhist(8,:), xhist(9,:))
title('Acceleration History','Interpreter','latex')
xlabel('x [km/sec\textsuperscript{2}]','Interpreter','latex')
ylabel('y [km/sec\textsuperscript{2}]','Interpreter','latex')
zlabel('z [km/sec\textsuperscript{2}]','Interpreter','latex')
grid on