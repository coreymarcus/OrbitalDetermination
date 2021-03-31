clear
close all
clc


%symbolic variables
syms r_x r_y r_z v_x v_y v_z 'real' %vehicle state ECI
syms Or_x Or_y Or_z 'real' 'real' %observer state ECI
syms w_earth 'real' %earth rotation rate
syms rel_x rel_y rel_z 'real' %relative ECI position
syms relv_x relv_y relv_z 'real' %relative ECI velocity

%get ECI velocity of the observer
Ov_x = -w_earth*Or_y;
Ov_y = w_earth*Or_x;
Ov_z = 0;

%position and velocity vectors
r = [r_x, r_y, r_z]';
v = [v_x, v_y, v_z]';
rO = [Or_x, Or_y, Or_z]';
vO = [Ov_x, Ov_y, Ov_z]';

%relative position and velocity
rel_pos = r - rO;
rel_vel = v - vO;
% rel_pos = [rel_x, rel_y, rel_z]';
% rel_vel = [relv_x, relv_y, relv_z]';

%range
rho = sqrt(rel_pos(1)^2 + rel_pos(2)^2 + rel_pos(3)^2);

%range rate
rho_dot = rel_pos'*rel_vel/rho;

%measurement
z = [rho, rho_dot]';

%jacobian
H = jacobian(z, [r; v]);

%make substitutions in H if able
H = subs(H,r_x-Or_x,rel_x);
H = subs(H,r_y-Or_y,rel_y);
H = subs(H,r_z-Or_z,rel_z);
H = subs(H,rel_vel(1),relv_x);
H = subs(H,rel_vel(2),relv_y); 
H = subs(H,rel_vel(3),relv_z); 