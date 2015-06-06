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
titleSize  = 25;
legendSize = 20;
axisSize   = 20;

subplot(rows,cols,1);
loglog(w,abs(cData)); hold on;
loglog(w,abs(G));
title('Magnitude','FontSize',titleSize);
legend('Orig','Levy','FontSize',legendSize);
xlabel('\omega','FontSize',axisSize);
ylabel('Mag (dB)','FontSize',axisSize);

subplot(rows,cols,2);
semilogx(w,rad2deg(phase(cData))); hold on;
semilogx(w,rad2deg(phase(G)));
title('Phase','FontSize',titleSize);
legend('Orig','Levy''FontSize',legendSize);
xlabel('\omega','FontSize',axisSize);
ylabel('\phi (deg)','FontSize',axisSize);

subplot(rows,cols,3);
loglog(w,abs(G)-abs(cData));
title('Magnitude Error','FontSize',titleSize);
xlabel('\omega','FontSize',axisSize);
ylabel('\Delta Mag (dB)','FontSize',axisSize);

subplot(rows,cols,4);
semilogx(w,rad2deg(phase(G))-rad2deg(phase(cData)));
title('Phase Error','FontSize',titleSize);
xlabel('\omega','FontSize',axisSize);
ylabel('\Delta \phi (deg)','FontSize',axisSize);

