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
x0 = [6980.3967323588
    1619.61802198332
    15.1399428739289
    -1.66690187359566
    7.2578409459164
    0.261907498000759
    0
    0
    0
    0
    0.020
    0];


%change units a bit
x0(4:9) = 1000*x0(4:9);

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display',...
    'iter','UseParallel',false,'StepTolerance',1E-8);

%function
fun = @(x) GetCost(x);

%optimize
xf = lsqnonlin(fun,x0,[],[],options);

%change units back
xf(4:9) = xf(4:9)/1000;

%write units
writematrix(xf,"../../data/optimized_x0.csv");
