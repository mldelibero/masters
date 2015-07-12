% scripts/params/lossTan.m
% Run this to generate a plot for the loss tangent of a capacitor

% Clear environment
clearvars;
close all;
format shorte;

% Add funs to path
addpath('../utilityFuns');
addpath('../');

theta = pi / 9; % 20 degrees
esr = [0,1];
jxc = esr ./ cos(theta);

[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.oneError, jxc, jxc, esr);
xlabel(ax,'ESR')
ylabel(ax,'jX_C');
h_title.String = 'Loss Tangent';
set(gca,'Xtick',[],'Ytick',[]);

%%
%figure;
xstart = 0.3 * max(esr);
r = 1.8*xstart;
x = linspace(xstart,r,10);
%xf = fliplr(x);
%xt = [x,xf];

y = [sqrt(r^2-x.^2)];%,-sqrt(r^2-xf.^2)];
hold on;
plot(x,y);

