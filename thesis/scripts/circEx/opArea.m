% Run this to generate a plot for the safe operating area for testing
clearvars;
close all;

w = 2*pi*logspace(-1,6,100);
C = 1*10^-6;

Zc = 1 ./ (1i .* w .* C);

[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.cData, Zc, Zc, w);
h_leg.Location = 'northeast';
h_title.String = 'Capacitor Impedance';

