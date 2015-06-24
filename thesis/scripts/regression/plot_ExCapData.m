% plot_ExCapData.m
% This script outputs a plot of the example data taken from Murata's Sim Surfing online tool.

% Clear environment
clearvars;
close all;
format shorte;

filename = './data/GRM31MR71H105KA88.txt';
[w, cData, rData, iData] = getData(filename);
myPlotType   = plotType.cData;
plotFit(myPlotType, cData, cData, w); % figures/modeling/levyIter.jpg

