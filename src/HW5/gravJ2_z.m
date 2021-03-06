function out1 = gravJ2_z(J2,Rearth,mu,x,y,z)
%GRAVJ2_Z
%    OUT1 = GRAVJ2_Z(J2,REARTH,MU,X,Y,Z)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    06-Apr-2021 10:01:00

t2 = x.^2;
t3 = y.^2;
t4 = z.^2;
out1 = J2.*Rearth.^2.*mu.*z.*1.0./(t2+t3+t4).^(7.0./2.0).*(t2.*3.0+t3.*3.0-t4.*2.0).*(-3.0./2.0);
