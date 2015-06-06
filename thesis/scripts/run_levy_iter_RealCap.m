% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
clearvars;
close all;
format shorte;

NumDeg = 7;
DenDeg = 7;
iterations = 100;
filename = './data/GRM31MR71H105KA88.txt';

[w, cData, rData, iData] = getData(filename);
G = regression_levy_iter(cData, w, iterations, NumDeg, DenDeg);

%%
% Plot
figure;
rows = 2;
cols = 2;

subplot(rows,cols,1);
semilogx(w,abs(cData)); hold on;
%semilogx(w,abs(G(1,1:size(G,2))))
semilogx(w,abs(G(iterations,1:size(G,2))))
title('Magnitude');
%legend('Orig','i1','i_n');
legend('Orig','i_n');

subplot(rows,cols,2);
semilogx(w,rad2deg(phase(cData))); hold on;
%semilogx(w,rad2deg(phase(G(1,1:size(G,2)))));
semilogx(w,rad2deg(phase(G(iterations,1:size(G,2)))));
title('Phase');
%legend('Orig','i1','i_n');
legend('Orig','i_n');

subplot(rows,cols,3);
%semilogx(w,abs(G(1,1:size(G,2)))-abs(cData)); hold on;
semilogx(w,abs(G(iterations,1:size(G,2)))-abs(cData));
%legend('i1','i_n');
legend('i_n');
title('Magnitude Error');

subplot(rows,cols,4);
%semilogx(w,rad2deg(phase(G(1,1:size(G,2))))-rad2deg(phase(cData))); hold on;
semilogx(w,rad2deg(phase(G(iterations,1:size(G,2))))-rad2deg(phase(cData)));
%legend('i1','i_n');
legend('i_n');
title('Phase Error');
