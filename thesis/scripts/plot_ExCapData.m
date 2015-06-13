% This script runs the example in Sanathanan's paper showing how to improve Levy's method by an iterative method.

% Clear environment
clearvars;
close all;
format shorte;

filename = './data/GRM31MR71H105KA88.txt';
[w, cData, rData, iData] = getData(filename);
myPlotType   = plotType.cData;
plotFit(myPlotType, cData, cData, w); % figures/modeling/levyIter.jpg

