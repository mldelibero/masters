function regression_highLevel()
    % Clear environment
    clearvars;
    close all;

    % Prepare the data
    if ~exist('filename') || ~exist('modelOrder')
        filename = './data/GRM31MR71H105KA88.txt';
        modelOrder  = 5;
    end
%    [w, cData, rData, iData] = getData(filename);
%    [T, w_norm] = calcChebyPoly(w);

    % Temporary data
     a = [1 2 3 2 1 4];
     b = [1 2 3 2 3];
  %   [cData,w] = freqs(b,a,64);
    w = logspace(0,6,1000)';
    [cData] = freqs(b,a,w);

    [numCoeff,denCoeff] = invfreqs(cData,w,4,5);
    G = calcOutput(fliplr(numCoeff), fliplr(denCoeff), w);
    abs_err = abs(cData) - abs(G);
    pha_err = phase(cData) - phase(G);

    figure;

    subplot(2,2,1);
    semilogx(w,abs(cData));
    hold on;
    semilogx(w,abs(G));
    title('Magnitude');

    subplot(2,2,2);
    semilogx(w,rad2deg(phase(cData)));
    hold on;
    semilogx(w,rad2deg(phase(G)));
    title('Phase');

    subplot(2,2,3);
    semilogx(w,abs_err);
    title('Magnitude Error');

    subplot(2,2,4);
    semilogx(w,rad2deg(pha_err));
    title('Phase Error');
end

function G = calcOutput(numCoeff, denCoeff, w)
    Num = numCoeff(1) * ones(length(w),1);
    Den = denCoeff(1) * ones(length(w),1);

    for h = 2:length(numCoeff);
        Num = Num + numCoeff(h) .* (1i * w).^(h-1);
    end

    for h = 2:length(denCoeff);
        Den = Den + denCoeff(h) .* (1i * w).^(h-1);
    end

    G = Num ./ Den;
end

