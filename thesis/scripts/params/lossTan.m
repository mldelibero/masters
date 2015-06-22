% Run this to generate a plot for the loss tangent of a capacitor
clearvars;
close all;

theta = pi / 9; % 20 degrees
esr = [0,1];
jxc = esr ./ cos(theta);

[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.oneError, jxc, jxc, esr);
xlabel(ax,'ESR')
ylabel(ax,'jX_C');
h_title.String = 'Loss Tangent';
set(gca,'Xtick',[],'Ytick',[]);

