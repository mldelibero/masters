% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
clearvars;
close all;
format shorte;

NumDeg = 3;
DenDeg = 3;
iterations = 0;
filename = './data/GRM31MR71H105KA88.txt';

[w, cData, rData, iData]      = getData(filename);
[G, numCoeffs, denCoeffs]     = regression_levy_iter(cData, w, iterations, NumDeg, DenDeg);
%%
% Plot
figure;
rows = 2;
cols = 2;

subplot(rows,cols,1);
loglog(w,abs(cData)); hold on;
loglog(w,abs(G));
title('Magnitude');
legend('Orig','Levy');

subplot(rows,cols,2);
semilogx(w,rad2deg(phase(cData))); hold on;
semilogx(w,rad2deg(phase(G)));
title('Phase');
legend('Orig','Levy');

subplot(rows,cols,3);
semilogx(w,abs(G)-abs(cData)); hold on;
title('Magnitude Error');

subplot(rows,cols,4);
semilogx(w,rad2deg(phase(G))-rad2deg(phase(cData))); hold on;
title('Phase Error');

