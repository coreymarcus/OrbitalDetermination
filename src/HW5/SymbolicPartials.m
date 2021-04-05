% Find partial derivatives of gravity potential equation using symbolic
% toolbox

clear
close all
clc

% constants
syms mu J2 J3 Rearth real
assume(mu > 0)
assume(J2 > 0)
assume(Rearth > 0)

%variables
syms x y z v_x v_y v_z Cd Csolar g real
R = [x y z]';
r = sqrt(R'*R);

%gravity potential
U_pm = mu/r;
U_J2 = -mu/r*J2*(Rearth/r)^2*(1.5*(z/r)^2 - 0.5);
U_J3 = mu*J3*(1/r^7)*0.5*z*(5*z^2 - 3*r^2);

%acceleration from gravity
gravPM = simplify(jacobian(U_pm,R));
gravJ2 = simplify(jacobian(U_J2,R));
gravJ3 = simplify(jacobian(U_J3,R));

%convert to files
matlabFunction(gravJ2(1),'File','gravJ2_x.m');
matlabFunction(gravJ2(2),'File','gravJ2_y.m');
matlabFunction(gravJ2(3),'File','gravJ2_z.m');

%find A
X = [x y z v_x v_y v_z]';
f = [v_x; v_y; v_z; gravPM'+gravJ2'+gravJ3'];
A = jacobian(f,X);
A_vec = reshape(A(4:6,1:3),9,1);
matlabFunction(A_vec,'File','createA.m');
