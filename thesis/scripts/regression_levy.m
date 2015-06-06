function G = regression_levy(cData, w, modelOrder)
    % Execute the regression analysis based upon Levy's method.
    [numCoeff, denCoeff] = calcCoeffs(w, cData, modelOrder);
    numCoeff
    denCoeff = [1;denCoeff]
    [G, Num, Den] = calcDataFromCoeffs(numCoeff, denCoeff, w);
end

function [TnumCoeff, TdenCoeff] = calcCoeffs(w, cData, modelOrder)
    % This calculates the polynomial coefficients from the wshev polynomial series

    % Prepare variables
    numSamples = size(w,1);
    matrixSize = modelOrder * 2 + 1;

    M = ones(matrixSize, matrixSize);
    N = ones(matrixSize,1);
    C = ones(matrixSize,1);

    sign  = 0;
    power = 0;

    format = 'Row:%i, Col:%i, Q:%i, P:%i, S:%i, F:%s\n';

    %% Calculate
    % Populate M
    for row = 1:matrixSize
        for col = 1:matrixSize
            % Find out which section of M you are in:
            if ((row <= matrixSize - modelOrder) && (col <= matrixSize - modelOrder))
                Q = 0;
                %Upper Left Quadrant -- lambdas
                if (xor(mod(row,2) == 0, mod(col,2) == 0) == 1) % If both are odd or both are even
                    M(row,col) = 0;
                    power = 0;
                    sign = 0;
                    F = '0';
                    %fprintf(format,row,col,Q,power,sign,F);
                else
                    power = row + col - 2;
                    if ((mod(col,3) == 0) && (mod(row,2) == 1))
                        sign = -1;
                    elseif ((mod(col,4) == 0) && (mod(row,2) == 0))
                        sign = -1;
                    else
                        sign = 1;
                    end
                    
                    M(row,col) = sign * lambda(w,power);
                    F = 'L';
                    %fprintf(format,row,col,Q,power,sign,F);
                end % if (xor(mod(row,2) == 0, mod(col,2) == 0) == 1) % If both are odd or both are even
            elseif (row <= matrixSize - modelOrder)
                % Upper Right quadrant -- S and T
                Q = 1;
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
                    M(row,col) = sign * T(w,imag(cData),power);
                    F = 'T';
                else
                    M(row,col) = sign * S(w,real(cData),power);
                    F = 'S';
                end

                %fprintf(format,row,col,Q,power,sign,F);
                
            elseif (col <= modelOrder+1)
                % Lower Left -- S and T
                Q = 2;
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
                    M(row,col) = sign * S(w,real(cData),power);
                    F = 'S';
                else
                    M(row,col) = sign * T(w,imag(cData),power);
                    F = 'T';
                end

                %fprintf(format,row,col,Q,power,sign,F);
            else
                % Lower Right -- U
                Q = 3;
                locrow = row - matrixSize + modelOrder;
                loccol = col - matrixSize + modelOrder;
                if (xor(mod(locrow,2) == 0, mod(loccol,2) == 0) == 1) % If both are odd or both are even
                    M(row,col) = 0;
                    power = 0;
                    sign = 0;
                    F = '0';
                    %fprintf(format,row,col,Q,power,sign,F);
                else
                    power = locrow + loccol;
                    if     ((mod(locrow,2) == 1) && (mod(loccol,3) == 0))
                        sign = -1;
                    elseif ((mod(locrow,2) == 0) && (mod(loccol,4) == 0))
                        sign = -1;
                    else
                        sign = 1;
                    end
                    
                    M(row,col) = sign * U(w,real(cData),imag(cData),power);
                    F = 'U';
                    %fprintf(format,row,col,Q,power,sign,F);
                end
            end
        end
    end
    % Populate C
    for (row = 1:matrixSize)
        if (row <= matrixSize - modelOrder)
            power = row - 1;
            if (mod(power,2) == 0)
                C(row) = S(w,real(cData),power);
            else
                C(row) = T(w,imag(cData),power);
            end
        else
            localRow = row - matrixSize + modelOrder;
            if (mod(localRow,2) == 1)
                C(row) = 0;
            else
            C(row) = U(w,real(cData),imag(cData),localRow);
            end
        end
    end
    
    % Solve For N
    %M(1:size(M,2),1)
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

