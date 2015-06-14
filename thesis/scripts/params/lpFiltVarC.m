% Run this to generate a plot for an RC low pass filter with a varying capacitance.
close all;

w = 2*pi*logspace(-1,6,100);

C1 = 1*10^-6;
C2 = 1*10^-5;
C3 = 1*10^-4;

R = 1000;
Zc = [1 ./ (1i*w*C1); 1 ./ (1i*w*C2); 1./ (1i*w*C3)];
Vout = Zc ./ (R + Zc);

plotFit(plotType.MULTCDATA, Vout, Vout, w)

