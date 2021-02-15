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
U_pm = mu/r;
U_J2 = -mu/r*J2*(Rearth/r)^2*(1.5*(z/r)^2 - 0.5);

%acceleration from gravity
gravPM = simplify(jacobian(U_pm,R));
gravJ2 = simplify(jacobian(U_J2,R));

%convert to files
matlabFunction(gravJ2(1),'File','gravJ2_x.m')
matlabFunction(gravJ2(2),'File','gravJ2_y.m')
matlabFunction(gravJ2(3),'File','gravJ2_z.m')

