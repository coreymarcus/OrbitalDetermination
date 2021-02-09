%This script does all the plotting of my data, because plotting in cpp is hard :/

%startup
clear
close all
clc

%savepath for images
savepath = '../../doc/HW1/figs/';

%load the data
xhist = csvread('../../data/xhist_HW1.csv');
t = csvread('../../data/thist_HW1.csv');
E = csvread('../../data/ehist_HW1.csv');

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

%calculate magnitudes
N = length(t);
posnorm = zeros(1,N);
velnorm = zeros(1,N);
accelnorm = zeros(1,N);
for ii = 1:N
    posnorm(ii) = norm(xhist(1:3,ii));
    velnorm(ii) = norm(xhist(4:6,ii));
    accelnorm(ii) = norm(xhist(7:9,ii));
end

figure
subplot(3,1,1)
plot(t, posnorm)
title('Magnitude of Position and Derivatives','Interpreter','latex')
ylabel('Radius [km]','Interpreter','latex')
grid on

subplot(3,1,2)
title('Magnitude of Position and Derivatives','Interpreter','latex')
plot(t, velnorm)
ylabel('Velocity [km/sec]','Interpreter','latex')
grid on

subplot(3,1,3)
title('Magnitude of Position and Derivatives','Interpreter','latex')
plot(t, accelnorm)
ylabel('Acceleration [km/sec\textsuperscript{2}]','Interpreter','latex')
grid on
xlabel('Time [sec]','Interpreter','latex')
saveas(gcf,[savepath 'posandderivs.pdf'])

%scatter angular momentum
h = xhist(10:12,:);
figure
scatter3(h(1,:),h(2,:),h(3,:),2,'filled')
title('Specific Angluar Momentum Vector','Interpreter','latex')
xlabel('x [km\textsuperscript{2}/sec]','Interpreter','latex')
ylabel('y [km\textsuperscript{2}/sec]','Interpreter','latex')
zlabel('z [km\textsuperscript{2}/sec]]','Interpreter','latex')
saveas(gcf,[savepath 'angmom.pdf'])

%plot energy
figure
subplot(2,1,1)
plot(t,E(1,:))
ylabel('KE [J/kg]','Interpreter','latex')
title('Specific Kinetic and Potential Energy','Interpreter','latex')

subplot(2,1,2)
plot(t,E(2,:))
ylabel('PE [J/kg]','Interpreter','latex')
xlabel('Time [sec]','Interpreter','latex')
saveas(gcf,[savepath 'energy.pdf'])

E0 = sum(E(:,1));

figure
plot(t,sum(E,1) - E0)
xlabel('Time [sec]', 'Interpreter', 'latex')
ylabel('$dE = E(t) - E(t_0)$', 'Interpreter', 'latex')
saveas(gcf,[savepath 'energychange.pdf'])
title('Change In Energy Over Time', 'Interpreter', 'latex')