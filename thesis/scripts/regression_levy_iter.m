function [A,M,N,C,G, numCoeffs, denCoeffs] = regression_levy_iter(cData, w, iterations, numbNumCoeffs, numbDenCoeffs)
    % Execute the regression analysis based upon Levy's method.
    % G(s) = (a0 + a1*s + a2*s^2 + ... + an*s^n) /
    %        (b0 + b1*s + b2*s^2 + ... + bm*s^m)

    % Solve
%    [TnumCoeff,TdenCoeff] = calcCoeffs(w_norm, T, cData, modelOrder);
    numCoeffs = ones(numbNumCoeffs,1);
    denCoeffs = ones(numbDenCoeffs,1); % Needs to be ones for the initial guess of W
    W         = ones(length(w)); % Weighting function

    %Calculate the Den for W from the initial guess of denCoeffs
    [Gtemp,Num,Den] = calcDataFromCoeffs(numCoeffs,denCoeffs,w);

    % Run Levy's original formula if not using iterations.
    if iterations == 0
        iterations = 1;
        Den = ones(length(Den));
    end

    for (iter = 1:iterations)
        W = 1 ./ abs(Den).^2;
        [A,M,N,C,numCoeffs, denCoeffs] = calcCoeffs(cData, w, W, numbNumCoeffs, numbDenCoeffs);
        [G(iter,1:length(w)),Num,Den] = calcDataFromCoeffs(numCoeffs,denCoeffs,w);
        iter
    end
end % function [G,numCoeffs,denCoeffs] = regression_levy_iter(cData, w, iterations, numbNumCoeffs, numbDenCoeffs)

function [A,M,N,C,numCoeffs, denCoeffs] = calcCoeffs(cData, w, W, numbNumCoeffs, numbDenCoeffs)
    % This calculates the polynomial coefficients from the wshev polynomial series

    % Prepare variables
    Mrows = numbNumCoeffs + numbDenCoeffs - 1;
    Mcols = Mrows;

    M = ones(Mrows,Mcols);
    N = ones(Mrows,1);
    C = ones(Mrows,1);
    A = -13*ones(Mcols*Mrows,8);

    sign  = 0;
    power = 0;

%    formatPrint = 'Row: %i, Col: %i, xRow: %i, xCol: %i, Q: %i, P: %i, S: %i\n';

    %% Calculate
    % Populate M
    for row = 1:Mrows
        for col = 1:Mcols
            % Find out which section of M you are in:
            if ((row <= numbNumCoeffs) && (col <= numbNumCoeffs))
                %Upper Left Quadrant -- lambdas
                xrow = row;
                xcol = col;
                Q    = 0;

                if (xor(mod(xrow,2) == 0, mod(xcol,2) == 0) == 1) % If both are odd or both are even
                    M(xrow,xcol) = 0;
                    power = 0;
                    sign = 0;
                    F = 0;
                   %fprintf(formatPrint,row,col,0);
                    A((row-1)*(Mcols-1)+row+col-1,:) = [row,col,xrow,col,Q,power,sign,F];

                else
                    power = xrow + xcol - 2;

                    if ((mod(xcol,3) == 0) && (mod(xrow,2) == 1))
                        sign = -1;
                    elseif ((mod(xcol,4) == 0) && (mod(xrow,2) == 0))
                        sign = -1;
                    else
                        sign =  1;
                    end
                        
                    M(xrow,xcol) = sign * lambda(W,w,power);
                    F = 1;
                end % else

                %fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);
                    A((row-1)*(Mcols-1)+row+col-1,:) = [row,col,xrow,col,Q,power,sign,F];

            elseif (row <= numbNumCoeffs)
                % Upper Right quadrant -- S and T
                xrow = row;
                xcol = col - numbNumCoeffs; 
                Q    = 1;

                power = xrow + xcol - 1;

                if     (mod(xcol,4) == 1) && (mod(xrow,2) == 1)
                    sign =  1;
                elseif (mod(xcol,4) == 1) && (mod(xrow,2) == 0)
                    sign = -1;
                elseif (mod(xcol,4) == 2)
                    sign =  1;
                elseif (mod(xcol,4) == 3) && (mod(xrow,2) == 1)
                    sign = -1;
                elseif (mod(xcol,4) == 3) && (mod(xrow,2) == 0)
                    sign = 1;
                elseif (mod(xcol,4) == 4)
                    sign = -1;
                end

                if  (mod(power,2) == 1)
                    M(row,col) = sign * T(W,w,imag(cData),power);
                    F = 3;
                else
                    M(row,col) = sign * S(W,w,real(cData),power);
                    F = 2;
                end

               %fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);
                    A((row-1)*(Mcols-1)+row+col-1,:) = [row,col,xrow,col,Q,power,sign,F];
                
            elseif ((row > numbNumCoeffs) && (col <= numbNumCoeffs))
                % Lower Left -- S and T
                xrow = row - numbNumCoeffs;
                xcol = col;
                Q    = 2;

                power = xrow + col - 1;

                if     (mod(xcol,4) == 1)
                    sign =  1;
                elseif (mod(xcol,4) == 2) && (mod(xrow,2) == 1)
                    sign = -1;
                elseif (mod(xcol,4) == 2) && (mod(xrow,2) == 0)
                    sign =  1;
                elseif (mod(xcol,4) == 3)
                    sign = -1;
                elseif (mod(xrow,2) == 1)
                    sign =  1;
                elseif (mod(xrow,2) == 0)
                    sign = -1;
                end
                
                if (mod(power,2) == 0)
                    M(row,col) = sign * S(W,w,real(cData),power);
                    F = 2;
                else
                    M(row,col) = sign * T(W,w,imag(cData),power);
                    F = 3;
                end
               %fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);
                    A((row-1)*(Mcols-1)+row+col-1,:) = [row,col,xrow,col,Q,power,sign,F];

            else
                % Lower Right -- U
                xrow = row - numbNumCoeffs;
                xcol = col - numbNumCoeffs;
                Q    = 3;

                if (xor(mod(xrow,2) == 0, mod(xcol,2) == 0) == 1) % If both are odd or both are even
                    M(row,col) = 0;
                    power = 0;
                    sign = 0;
                    F = 0;
                    A((row-1)*(Mcols-1)+row+col-1,:) = [row,col,xrow,col,Q,power,sign,F];
                else
                    power = xrow + xcol;

                    if     (mod(xcol,4) == 1)
                        sign = 1;
                    elseif (mod(xcol,4) == 2) 
                        sign = 1;
                    else
                        sign = -1;
                    end

                    M(row,col) = sign * U(W,w,real(cData),imag(cData),power);
                    F = 4;
                end % else -- check for zeros
                 %fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);
                    A((row-1)*(Mcols-1)+row+col-1,:) = [row,col,xrow,col,Q,power,sign,F];
            end % else -- searching through quadrants
        end % for row = 1:Mrows
    end % for col = 1:Mcols

    % Populate C
   %fprintf('Calc C\n');
   %fprintf('------\n');
    cFormat = 'row:%i, xrow:%i, P:%i, F:%s\n';
    for (row = 1:length(C))
        if (row <= numbNumCoeffs)
            xrow = row;
            power = xrow - 1;
            if (mod(power,2) == 0)
                C(xrow) = S(W,w,real(cData),power);
               %fprintf(cFormat,row,row,power,2);
            else
                C(xrow) = T(W,w,imag(cData),power);
               %fprintf(cFormat,row,xrow,power,3);
            end
        else
            xrow = row - numbNumCoeffs;
            power = xrow;
            if (mod(xrow,2) == 1)
                C(row) = 0;
               %fprintf(cFormat,row,xrow,power,0);
            else
                C(row) = U(W,w,real(cData),imag(cData),power);
               %fprintf(cFormat,row,xrow,power,4);
            end
        end
    end
    
    % Solve For N
    N = inv(M) * C;

    numCoeffs = N(1:numbNumCoeffs);
    denCoeffs = [1;N(numbNumCoeffs+1:length(N))]; % Need to add 1 as coeff_0
 
    %Sub Functions
    function [Lh] = lambda(W,wk,h)
        Lh = 0;
        for k = 1:1:length(wk)
            Lh = Lh + wk(k)^h*W(k);
        end
    end

    function [S] = S(W,wk,Rk,h)
        S = 0;
        for k = 1:1:length(wk)
            S = S + wk(k)^h * Rk(k)*W(k);
        end
    end

    function [T] = T(W,wk,Ik,h)
        T = 0;
        for k = 1:1:length(wk)
            T = T + wk(k)^h * Ik(k)*W(k);
        end
    end

    function [U] = U(W,wk,Rk,Ik,h)
        U = 0;
        for k = 1:1:length(wk)
            U = U + wk(k)^h * (Rk(k)^2 + Ik(k)^2)*W(k);
        end
    end
end % function [b,a] = calcCoeffs(w_norm, T, cData)

