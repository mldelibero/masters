% run_basicLse.m
% This script runs the basic LSE for a line

% Clear environment
clearvars;
close all;
format shorte;

%% Generate data
n = 5;
err = rand(1,n)*1.3;
x = linspace(1,n,n);
b = rand(1)*ones(1,n);
y = x + err + b;

%% LSE
avg_y   = sum(y);
avg_x   = sum(x);
avg_xy  = sum(y .* x);
avg_xsq = sum(x .^ 2);

a_1 = (avg_xy - avg_x *avg_y) / (avg_xsq - avg_x^2)
a_0 = avg_y - a_1 * avg_x

y_lse = a_0*ones(1,length(x)) + a_1 .* x;

%% Plot
plotFit(plotType.DIFF_PLOT,y,y_lse,x);

