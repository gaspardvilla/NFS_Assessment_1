%% Clean
close all
clearvars
clc


%% Initialization of the general parameters

% Initailization of the bounds
lower = [-inf -inf -inf -1 0];
upper = [inf inf inf 1 inf];

% Initialization of the function
rmse_fct = @(x) RMSE(x);

% Initial values for [theta kappa sigma rho V0]
x0 = [0.04 1.5 0.3 -0.6 0.0441];


[x_min, f_min, exitflag, output] = fminsearchcon(rmse_fct, x0, lower, upper);


