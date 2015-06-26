% scripts/appendix/cheby/powerSeries.m
% This script plots a normalized power series
% It duplicates Fig1 of Reference: \cite{beyene_uwave}

% Clear environment
clearvars;
close all;
format shorte;

% Add funs to path
addpath('../../utilityFuns');
addpath('../../');

x = linspace(0,1,100);
order = 25;
w = ones(order,length(x));

x = x ./ max(x); % Normalize x
w(2,:) = x;
%w(2,:) = x ./ max(x);

for index = 3:order
    w(index,:) = x.^(index-1);
%    w(index,:) = w(index,:) ./ max(w(index,:));
end

[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.MULTDATA, w, w, x); % figures/appendix/cheby/cheby.jpg
h_title.String = 'Power Series';
h_leg = legend(p(1:5),'w^0','w^1','w^2','w^3','w^4');
xlabel(ax,'x');
ylabel(ax,'Power Series');

