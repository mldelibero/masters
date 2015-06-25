% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
clearvars;
close all;
format shorte;

NumDeg = 3;
DenDeg = 3;
iterations = 1;
filename = './data/GRM31MR71H105KA88.txt';
modelType = modelTypes.NO_MODEL;

[w, cData, rData, iData]      = getData(filename);
initDen = getInitGuess(w,modelType);
[G, numCoeffs, denCoeffs]     = regression_levy_iter(cData, w, iterations, NumDeg, DenDeg, initDen);

%% Plot
plotcDiff = plotType.cVectorsDiff;
plotFit(plotcDiff, cData, G, w); % figures/regression/levy.jpg

