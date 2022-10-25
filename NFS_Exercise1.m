close all
clear all
clc


%% Picard's method

% Initialization of the constants
Ta = 20;
Tb = 40;
k = 400;

% Function S 
S = @(T) 4 - (5 * (T^3));
S_approx = @(T) 0;
