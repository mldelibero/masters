% This script runs the iterative regression analysis to fit a RSL capacitor model

% Clear environment
clearvars;
close all;
format shorte;

NumDeg = 3;
DenDeg = 1;
iterations = 100;
filename = './data/GRM31MR71H105KA88.txt';

[w, cData, rData, iData] = getData(filename);
[G, numCoeffs, denCoeffs] = regression_levy_iter(cData, w, iterations, NumDeg, DenDeg);

%% Error minimization
E = sumError(cData,G);
%%
Emag = E(:,1);
Epha = E(:,2);
EmagNorm = Emag ./ max(Emag);
EphaNorm = Epha ./ max(Epha);
E2 = EmagNorm + EphaNorm;
n = find(E2==min(E2),1);
fprintf('Error Minimized at Iteration: %i\n',n);
denCoeffs = [0;1];
[G2,Num,Den] = calcDataFromCoeffs(numCoeffs,denCoeffs,w);
a_0 = numCoeffs(1);
numCoeffs = numCoeffs ./ a_0;
denCoeffs = denCoeffs ./ a_0;
[G3,Num,Den] = calcDataFromCoeffs(numCoeffs,denCoeffs,w);

%% Only take what we want
% numbNumCoeffs = 3;
% numbDenCoeffs = 1;
% 
% numCoeffs2 = numCoeffs(1:numbNumCoeffs);
% denCoeffs2 = denCoeffs(1:numbDenCoeffs);
% 
% [G2,Num2,Den2] = calcDataFromCoeffs(numCoeffs2,denCoeffs2,w);
% % 
%% Plot 
plotcDiff = plotType.cVectorsDiff;
plotErrs  = plotType.twoErrors;
plotErr   = plotType.oneError;

plotFit(plotcDiff, cData, G(n,:), w); % figures/modeling/levyIter.jpg
plotFit(plotcDiff, cData, G2    , w); % figures/modeling/levyIter.jpg
plotFit(plotcDiff, cData, G3    , w); % figures/modeling/levyIter.jpg
%plotFit(plotErrs , Emag,  Epha     ); % figures/modeling/levyIter_Err1.jpg
%plotFit(plotErr  , E2              ); % figures/modeling/levyIter_Err2.jpg

