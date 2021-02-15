% Find partial derivatives of gravity potential equation using symbolic
% toolbox

clear
close all
clc

% constants
syms mu J2 Rearth real
assume(mu > 0)
assume(J2 > 0)
assume(Rearth > 0)

%variables
syms x y z real
R = [x y z]';
r = sqrt(R'*R);

%gravity potential
U = mu/r*(1 - J2*(Rearth/r)^2*(1.5*sin(z/r)^2 - 0.5));

%acceleration from gravity
grav = jacobian(U,R);

