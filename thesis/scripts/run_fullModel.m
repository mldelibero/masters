% This script runs the iterative regression analysis to fit a RSL capacitor model

% Clear environment
clearvars;
close all;
format shorte;

modelType = modelTypes.FULL_MODEL;
NumDeg = 4;
DenDeg = 3;
iterations = 1000;
filename = './data/GRM31MR71H105KA88.txt';


[w, cData, rData, iData] = getData(filename);
%[initDen] = getInitGuess(w,modelTypes.NO_MODEL);   % Uncomment to run the bad  output version
[initDen] = getInitGuess(w,modelTypes.FULL_MODEL); % Uncomment to run the good output version
[G, numCoeffs, denCoeffs, E, minIndex] = regression_levy_iter(cData, w, iterations, NumDeg, DenDeg, initDen);

%% Error minimization
Emag = E(:,1);
Epha = E(:,2);
EmagNorm = Emag ./ max(Emag);
EphaNorm = Epha ./ max(Epha);
E2 = EmagNorm + EphaNorm;
n = find(E2==min(E2),1);
fprintf('Error Minimized at Iteration: %i\n',n);

%% Calculate Circuit Parameters
%[LE, RE, RL, C, CD, RD] = calcCircParams(numCoeffs, denCoeffs);
params = calcCircParams(numCoeffs, denCoeffs, modelType);


%% Plot 
plotcDiff = plotType.cVectorsDiff;
plotErrs  = plotType.twoErrors;
plotErr   = plotType.oneError;

plotFit(plotcDiff, cData, G(n,:), w); % figures/modeling/levyIter.jpg
plotFit(plotErrs , Emag,  Epha     ); % figures/modeling/levyIter_Err1.jpg
plotFit(plotErr  , E2              ); % figures/modeling/levyIter_Err2.jpg

%% Plot Circuit Parameters
rows = 1;
cols = 2;
figure;
subplot(rows,cols,1); plot(params(1:10,1)); title('LE');
subplot(rows,cols,2); plot(params(1:10,2)); title('RE');

