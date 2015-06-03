function N = regression_levy(filename,modelOrder)
    % Execute the regression analysis based upon Levy's method.

    % Clear environment
    clearvars;
    close all;
    format shorte;

    % Prepare the data
    if ~exist('filename') || ~exist('modelOrder')
        filename = './data/GRM31MR71H105KA88.txt';
        modelOrder  = 5;
    end

    [w, cData, rData, iData] = getData(filename);
    [T, w_norm] = calcChebyPoly(w);

    % Solve
    [TnumCoeff,TdenCoeff] = calcCoeffs(w_norm, T, cData, modelOrder);
    [numCoeff] = clenshaw(TnumCoeff);
    [denCoeff] = clenshaw(TdenCoeff);

    % Plot
    plotOrigData(cData,w);
    plotAns(numCoeff,denCoeff,w);
end

function [TnumCoeff, TdenCoeff] = calcCoeffs(w_norm, cheby, cData, modelOrder)
    % This calculates the polynomial coefficients from the Chebyshev polynomial series

    % Prepare variables
    numSamples = size(cheby,1);
    matrixSize = modelOrder * 2 + 1;

    M = ones(matrixSize, matrixSize);
    N = ones(matrixSize,1);
    C = ones(matrixSize,1);

    sign  = 0;
    power = 0;

    %% Calculate
    % Populate M
    for row = 1:matrixSize
        for col = 1:matrixSize
            
            % Find out which section of M you are in:
            if ((row <= matrixSize - modelOrder) && (col <= matrixSize - modelOrder))
                %Upper Left Quadrant -- lambdas
                if (xor(mod(row,2) == 0, mod(col,2) == 0) == 1) % If both are odd or both are even
                    M(row,col) = 0;
                else
                    power = row + col - 2;
                    if ((mod(col,3) == 0) && (mod(row,2) == 1))
                        sign = -1;
                    elseif ((mod(col,4) == 0) && (mod(row,2) == 0))
                        sign = -1;
                    else
                        sign = 1;
                    end
                    
                    M(row,col) = sign * lambda(cheby,power);
            end
            elseif (row <= matrixSize - modelOrder)
                % Upper Right quadrant -- S and T
                power = col - matrixSize + modelOrder + row - 1;
                if (col == matrixSize - modelOrder + 1)
                    if (mod(power,2) == 0)
                        sign = -1;
                    else
                        sign = 1;
                    end
                elseif (col == matrixSize - modelOrder + 2)
                    sign = 1;
                elseif ((mod(row,2) == 1) && (mod(col - matrixSize + modelOrder,4) == 0))
                    sign = -1;
                elseif ((mod(row,2) == 0) && (mod(col - matrixSize + modelOrder,4) == 0))
                    sign = -1;
                else
                    sign = 1;
                end
                
                if  (mod(power,2) == 1)
                    M(row,col) = sign * T(cheby,imag(cData),power);
                else
                    M(row,col) = sign * S(cheby,real(cData),power);
                end
                
            elseif (col <= modelOrder+1)
                % Lower Left -- S and T
                power = col + row - matrixSize + modelOrder - 1;
                if (mod(col,3) == 0)
                    sign = -1;
                elseif (mod(col,3) == 1)
                    sign = 1;
                elseif ((mod(col,4) == 2) && mod(power,2) == 0)
                    sign = -1;
                elseif ((mod(col,4) == 0) && mod(power,2) == 1)
                    sign = -1;
                else
                    sign = 1;
                end
                
                if (mod(power,2) == 0)
                    M(row,col) = sign * S(cheby,real(cData),power);
                else
                    M(row,col) = sign * T(cheby,imag(cData),power);
                end
            else
                % Lower Right -- U
                locrow = row - matrixSize + modelOrder;
                loccol = col - matrixSize + modelOrder;
                if (xor(mod(locrow,2) == 0, mod(loccol,2) == 0) == 1) % If both are odd or both are even
                    M(row,col) = 0;
                else
                    power = locrow + loccol;
                    if     ((mod(locrow,2) == 1) && (mod(loccol,3) == 0))
                        sign = -1;
                    elseif ((mod(locrow,2) == 0) && (mod(loccol,4) == 0))
                        sign = -1;
                    else
                        sign = 1;
                    end
                    
                    M(row,col) = sign * U(cheby,real(cData),imag(cData),power);
                end
            end
        end
    end
    % Populate C
    for (row = 1:matrixSize)
        if (row <= matrixSize - modelOrder)
            power = row - 1;
            if (mod(power,2) == 0)
                C(row) = S(cheby,real(cData),power);
            else
                C(row) = T(cheby,imag(cData),power);
            end
        else
            localRow = row - matrixSize + modelOrder;
            if (mod(localRow,2) == 1)
                C(row) = 0;
            else
            C(row) = U(cheby,real(cData),imag(cData),localRow);
            end
        end
    end
    
    % Solve For N
    N = inv(M) * C;
    TnumCoeff = N(1:(matrixSize+1)/2);
    TdenCoeff = N(  (matrixSize+1)/2+1:matrixSize);
 
    %Sub Functions
    function [Lh] = lambda(wk,h)
        Lh = 0;
        for k = 1:1:length(wk)
            Lh = Lh + wk(k)^h;
        end
    end

    function [S] = S(wk,Rk,h)
        S = 0;
        for k = 1:1:length(wk)
            S = S + wk(k)^h * Rk(k);
        end
    end

    function [T] = T(wk,Ik,h)
        T = 0;
        for k = 1:1:length(wk)
            T = T + wk(k)^h * Ik(k);
        end
    end

    function [U] = U(wk,Rk,Ik,h)
        U = 0;
        for k = 1:1:length(wk)
            U = U + wk(k)^h * (Rk(k)^2 + Ik(k)^2);
        end
    end
end % function [b,a] = calcCoeffs(w_norm, T, cData)

function [coeff] = clenshaw(chebyCoeff)
    % This function transforms chebyshev coefficients into power series coefficients
    coeff = chebyCoeff;

    n  = length(coeff);
    y2 = zeros(n);
    y1 = zeros(n);
    y0 = zeros(n);
    
    %yk = 2*w*y1 - y2 + chebyCoeff[k)

    for k=1:n
        for index = n:-1:2
            %Left shift to multiply by w
            y0(index) = y1(index - 1);
        end
        yo = y0 .* 2 - y2;
        y0(1) = y0(1) + chebyCoeff(k);
    end
end % function [coeff] = clenshaw(chebyCoeff)

function plotOrigData(cData,w)
    % Plot
    figure;
    rows = 1;

    cols = 2;

    subplot(rows,cols,1);
    semilogx(w,abs(cData));
    title('Magnitude');

    subplot(rows,cols,2);
    semilogx(w,rad2deg(phase(cData)));
    title('Phase');
    end

function plotAns(numCoeff,denCoeff,w)
    % Calculate output from coefficients 
    Num = numCoeff(1) * ones(size(w,1),1);
    Den = ones(size(w,1),1);

    for h = 1:size(denCoeff);
        Num = Num + numCoeff(h+1) .* (1i * w).^h;
        Den = Den + denCoeff(h)   .* (1i * w).^h;
    end

    G = Num ./ Den;

    % Plot
    figure;
    rows = 1;
    cols = 2;

    subplot(rows,cols,1);
    semilogx(w,abs(G));
    title('Magnitude');

    subplot(rows,cols,2);
    semilogx(w,rad2deg(phase(G)));
    title('Phase');
end

