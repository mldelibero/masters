% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
%clearvars;
%close all;
format shorte;

NumDeg = 7;
DenDeg = 7;
iterations = 0;
filename = './data/GRM31MR71H105KA88.txt';

[w, cData, rData, iData] = getData(filename);
%[w_norm] = calc_wNorm(w) + 1; % Push range to [0,2]
w_max  = max(w);
w_norm = w ./ w_max;
[G, numCoeffs, denCoeffs] = regression_levy_iter(cData, w_norm, iterations, NumDeg, DenDeg);

for(coeff = 2:numel(numCoeffs))
    numCoeffs(coeff) = numCoeffs(coeff) / w_max^(coeff-1);
end
for(coeff = 2:numel(denCoeffs))
    denCoeffs(coeff) = denCoeffs(coeff) / w_max^(coeff-1);
end
[G, Num, Den] = calcDataFromCoeffs(numCoeffs,denCoeffs,w);

%%
% Plot
figure;
rows = 2;
cols = 2;

n1 = 1;
n2 = 1;
n3 = 1;
l1 = sprintf('iter:%i',n1);
l2 = sprintf('iter:%i',n2);
l3 = sprintf('iter:%i',n3);

subplot(rows,cols,1);
semilogx(w,abs(cData)); hold on;
semilogx(w,abs(G(n1,1:size(G,2))));
semilogx(w,abs(G(n2,1:size(G,2))));
semilogx(w,abs(G(n3,1:size(G,2))));
title('Magnitude');
legend('Orig',l1,l2,l3);

subplot(rows,cols,2);
semilogx(w,rad2deg(phase(cData))); hold on;
semilogx(w,rad2deg(phase(G(n1,1:size(G,2)))));
semilogx(w,rad2deg(phase(G(n2,1:size(G,2)))));
semilogx(w,rad2deg(phase(G(n3,1:size(G,2)))));
title('Phase');
legend('Orig',l1,l2,l3);

subplot(rows,cols,3);
semilogx(w,abs(G(n1,1:size(G,2)))-abs(cData)); hold on;
semilogx(w,abs(G(n2,1:size(G,2)))-abs(cData));
semilogx(w,abs(G(n3,1:size(G,2)))-abs(cData));
legend(l1,l2,l3);
title('Magnitude Error');

subplot(rows,cols,4);
semilogx(w,rad2deg(phase(G(n1,1:size(G,2))))-rad2deg(phase(cData))); hold on;
semilogx(w,rad2deg(phase(G(n2,1:size(G,2))))-rad2deg(phase(cData)));
semilogx(w,rad2deg(phase(G(n3,1:size(G,2))))-rad2deg(phase(cData)));
legend(l1,l2,l3);
title('Phase Error');

