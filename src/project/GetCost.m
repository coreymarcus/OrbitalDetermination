function [J] = GetCost(x)
%GetCost Matlab Wrapper for cost function

%fix units
% x(4:6) = x(4:6)/1000;
% x(7:12) = x(7:12)/1000000;

%item to write
xwrite = zeros(13,1);
% xwrite(1:6) = x(1:6);
% xwrite(13) = x(7);
xwrite(1:6) = [6978.83947078333
    1617.08566078009
    19.5045324160835
    -1.66346314624123
    7.26036443567597
    0.270402425183416];
xwrite(13) = x;

%write x to csv
writematrix(xwrite,"../../data/xeval_project.csv");

%change directory
cd("../../build");

%evaluate
system("./project_cost");

%return
cd("../src/project");

%read cost
J = csvread("../../data/evalcost_project.csv");
end

