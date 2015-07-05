% scripts/regression/run_levy_iter.m
% This script runs the example in Sanathanan's paper showing how 
% to improve Levy's method by an iterative method.

% Clear environment
clearvars;
close all;
format shorte;

% Add funs to path
addpath('../utilityFuns');
addpath('../');

NumDeg = 3;
DenDeg = 3;
iterations = 100;
filename = '../data/fullModel.txt';

[w, cData, rData, iData] = getData(filename);
initDen = getInitGuess(w,modelTypes.FULL_MODEL);
[G, numCoeffs, denCoeffs, E, minIndex] = ...
regression_levy_iter(cData, w, iterations, NumDeg, DenDeg, initDen);

%% Error minimization
Emag = E(:,1);
Epha = E(:,2);
EmagNorm = Emag ./ max(Emag);
EphaNorm = Epha ./ max(Epha);
E2 = EmagNorm + EphaNorm;
n = find(E2==min(E2),1);

%% Relative Error
abs_magErr = abs(G(minIndex,:)) - abs(cData);
abs_phaErr = phase(G(minIndex,:)) - phase(cData);

rel_magErr = abs(abs_magErr ./ abs(cData))*100;
rel_phaErr = abs(abs_phaErr ./ phase(cData))*100;

%% Plot 
plotcDiff = plotType.cVectorsDiff;
plotErrs  = plotType.twoErrors;
plotErr   = plotType.oneError;
plotRelErr= plotType.MULTPLOT;

plotFit(plotcDiff , cData, G(n,:), w); % figures/regression/levyIter.jpg
plotFit(plotErrs  , Emag , Epha); % figures/regression/levyIter_Err1.jpg
plotFit(plotErr   , E2         ); % figures/regression/levyIter_Err2.jpg
[h_fig, h_leg, h_title, ax, p] = plotFit(plotRelErr, rel_magErr, rel_phaErr,w);
h_title(1).String = 'Relative Magnitude Error';
h_title(2).String = 'Relative Phase Error';
xlabel(ax(1),'w (rad)');
xlabel(ax(2),'w (rad)');
ylabel(ax(1,1), 'Mag (%)');
ylabel(ax(1,2), '\phi (%)');


