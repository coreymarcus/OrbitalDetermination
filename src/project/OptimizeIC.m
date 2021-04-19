clear
close all
clc


%initial guess
% x0 = [6984.45711518852;
%     1612.2547582643;
%     13.0925904314402;
%     -1.67667852227336;
%     7.26143715396544;
%     0.259889857225218
%     0
%     0
%     0
%     0
%     0.020
%     0];
% x0 = [6978.83947078333
%     1617.08566078009
%     19.5045324160835
%     -1.66346314624123
%     7.26036443567597
%     0.270402425183416
%     0
%     0
%     0
%     0
%     0
%     0
%     1.88];
% x0 = [6978.83947078333
%     1617.08566078009
%     19.5045324160835
%     -1.66346314624123
%     7.26036443567597
%     0.270402425183416
%     1.88];
x0 = [1.88];


%change units a bit
% x0(4:6) = 1000*x0(4:6);
% x0(7:12) = 1000000*x0(7:12);

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display',...
    'iter-detailed','UseParallel',false,'StepTolerance',1E-8,...
    'MaxFunctionEvaluations',1E8,'MaxIterations',1E4);

%function
fun = @(x) GetCost(x);

%optimize
xf = lsqnonlin(fun,x0,[],[],options);

%change units back
% xf(4:6) = xf(4:6)/1000;
% xf(7:12) = xf(7:12)/1000000;

%write units
writematrix(xf,"../../data/optimized_x0.csv");
