% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
clearvars;
%close all;
format shorte;

NumDeg = 6;
DenDeg = 7;

%numCoeff = [1.2768, 1.2803, 7.8236*10^-1, 7.9196*10^-2, 8.0901*10^-3, 2.8952*10^-4, 2.0144*10^-5];
%denCoeff = [1.0000, 2.5313, 4.2704*10^-1, 5.4648*10^-2, 4.5377*10^-3, 1.9841*10^-4, 1.1451*10^-5, 1.0522*10^-7];
numCoeff = [0.99936, 1.0086, -0.000015983];
denCoeff = [1.00000, 0.10097, 0.010031];
%w = linspace(.3,40,1000);
w = logspace(-1,2,100);

H = calcDataFromCoeffs(numCoeff,denCoeff,w);
G = regression_levy_iter(H, w, length(numCoeff), length(denCoeff));

% Plot
%figure;
rows = 2;
cols = 2;

subplot(rows,cols,1);
semilogx(w,abs(H)); hold on;
semilogx(w,abs(G));
title('Magnitude');
legend('Orig','Calc');

subplot(rows,cols,2);
semilogx(w,rad2deg(phase(H))); hold on;
semilogx(w,rad2deg(phase(G)));
title('Phase');
legend('Orig','Calc');

subplot(rows,cols,3);
semilogx(w,abs(G) - abs(H));
title('Magnitude Error');

subplot(rows,cols,4);
semilogx(w,rad2deg(phase(G)) - rad2deg(phase(H)));
title('Phase Error');

