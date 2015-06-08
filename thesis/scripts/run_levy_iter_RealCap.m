% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
clearvars;
close all;
format shorte;

NumDeg = 7;
DenDeg = 7;
iterations = 1000;
filename = './data/GRM31MR71H105KA88.txt';

[w, cData, rData, iData] = getData(filename);
[G, numCoeffs, denCoeffs] = regression_levy_iter(cData, w, iterations, NumDeg, DenDeg);

%% Error minimization
E = sumError(cData,G);
%%
Emag = E(:,1);
Epha = E(:,2);
EmagNorm = Emag ./ max(Emag);
EphaNorm = Epha ./ max(Epha);
E2 = EmagNorm + EphaNorm;
n = find(E2==min(E2),1);
fprintf('Error Minimized at Iteration: %i\n',n);

%% Plot 
plotcDiff = plotType.cVectorsDiff;
plotErrs  = plotType.twoErrors;
plotErr   = plotType.oneError;

plotFit(plotcDiff, cData, G(n,:), w); % figures/modeling/levyIter.jpg
plotFit(plotErrs , Emag,  Epha     ); % figures/modeling/levyIter_Err1.jpg
plotFit(plotErr  , E2              ); % figures/modeling/levyIter_Err2.jpg

