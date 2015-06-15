% Run this to generate a plot for the impedance of an LC capacitor. 
clearvars;
close all;

w = 2*pi*logspace(-1,6,1000);
C = 1*10^-6;
L = 1*10^-3;

Z = 1 ./ (1i .* w .* C) + 1i.*w.*L;

[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.cData, Z, Z, w);
h_leg.Location = 'northeast';
h_title.String = '(ESL) Capacitor Impedance';

