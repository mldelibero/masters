% scripts/circEx/opArea.m
% This script plots the Safe/Usable operating area for the circuit. 

% Clear environment
clearvars;
close all;
format shorte;

% Add funs to path
addpath('../utilityFuns');
addpath('../');

%% Set Variables
R = [1.1*10^3;1.3*10^5;1*10^7];
Tau = 0.1 * ones(length(R),1);
C = Tau ./ R .* 10^6;

P = [3;5;3];
R2 = [R(1);R(1);R;R(3);R(3);flipud(R)];
V = [0;sqrt(R2(1)*P(1));sqrt(R(1)*P(1));500;500;500;0;0;0;0];

%% Plot
%semilogx(R2,V); hold on;
[h_fig, h_leg, h_title, ax, p] = plotFit(plotType.OPAREA , V, C, R2); % figures/circEx/opArea.jpg

