function [J] = GetCost(x)
%GetCost Matlab Wrapper for cost function

%fix units
x(4:9) = x(4:9)/1000;

%write x to csv
writematrix(x,"../../data/xeval_project.csv");

%change directory
cd("../../build");

%evaluate
system("./project_cost");

%return
cd("../src/project");

%read cost
J = csvread("../../data/evalcost_project.csv");
end

