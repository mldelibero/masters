% scripts/appendix/cheby/cheby.m
% This script plots several of the first Chebyshev ploynomials

% Clear environment
clearvars;
close all;
format shorte;

% Add funs to path
addpath('../../utilityFuns');
addpath('../../');

x = linspace(-1,1,1000);
order = 6;
T = ones(order,length(x));

T(2,:) = T(2,:).* x;

for index = 3:order
    T(index,:) = 2 .* x .* T(index-1,:) - T(index-2,:);
end

[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.MULTDATA, T, T, x); % figures/appendix/cheby/cheby.jpg
h_title.String = 'Chebyshev Polynomials';
h_leg = legend(p,'T_0(x)','T_1(x)','T_2(x)','T_3(x)','T_4(x)','T_5(x)');
ylabel(ax,'T');

