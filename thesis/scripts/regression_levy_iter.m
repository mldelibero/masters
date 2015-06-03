function G = regression_levy_iter(cData, w, numbNumCoeffs, numbDenCoeffs)
    % Execute the regression analysis based upon Levy's method.
    % G(s) = (a0 + a1*s + a2*s^2 + ... + an*s^n) /
    %        (b0 + b1*s + b2*s^2 + ... + bm*s^m)

    % Solve
%    [TnumCoeff,TdenCoeff] = calcCoeffs(w_norm, T, cData, modelOrder);
    numCoeffs = ones(length(numbNumCoeffs));
    denCoeffs = ones(length(numbDenCoeffs));
    W         = zeros(length(w)); % Weighting function
    iters = 1;

    for (iter = 1:iters)
        [G,Num,Den] = calcDataFromCoeffs(numCoeffs,denCoeffs,w);
        W = 1 ./ abs(Den).^2;
        [numCoeffs, denCoeffs] = calcCoeffs(cData, w, W, numbNumCoeffs, numbDenCoeffs);
    end

    G = calcDataFromCoeffs(numCoeffs,denCoeffs,w);
end

function [numCoeffs, denCoeffs] = calcCoeffs(cData, w, W, numbNumCoeffs, numbDenCoeffs)
    % This calculates the polynomial coefficients from the wshev polynomial series

    % Prepare variables
    Mrows = numbNumCoeffs + numbDenCoeffs - 1;
    Mcols = Mrows;

    M = ones(Mrows,Mcols);
    N = ones(Mrows,1);
    C = ones(Mrows,1);

    sign  = 0;
    power = 0;

%%    formatPrint = 'Row: %i, Col: %i, xRow: %i, xCol: %i, Q: %i, P: %i, S: %i\n';
    
%    formatPrint = 'Row: %i, Col: %i, F: %s\n';

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
%                    fprintf(formatPrint,row,col,'0');
                else
%                    fprintf(formatPrint,row,col,'L');
                    power = xrow + xcol - 2;

                    if ((mod(xcol,3) == 0) && (mod(xrow,2) == 1))
                        sign = -1;
                    elseif ((mod(xcol,4) == 0) && (mod(xrow,2) == 0))
                        sign = -1;
                    else
                        sign =  1;
                    end
                        
                    M(xrow,xcol) = sign * lambda(w,power);
                end % else

%%                fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);

            elseif (row <= numbNumCoeffs)
                % Upper Right quadrant -- S and T
                xrow = row;
                xcol = col - numbNumCoeffs; 
                Q    = 1;

                power = row + xcol - 1;

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
                    M(row,col) = sign * T(w,imag(cData),power);
%                    fprintf(formatPrint,row,col,'T');
                else
                    M(row,col) = sign * S(w,real(cData),power);
%                    fprintf(formatPrint,row,col,'S');
                end

%%                fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);
                
            elseif ((row > numbNumCoeffs) && (col <= numbNumCoeffs))
                % Lower Left -- S and T
                xrow = row - numbNumCoeffs;
                xcol = col;
                Q    = 3;

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
                    M(row,col) = sign * S(w,real(cData),power);
%                    fprintf(formatPrint,row,col,'S');
                else
                    M(row,col) = sign * T(w,imag(cData),power);
%                    fprintf(formatPrint,row,col,'T');
                end
%%                fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);

            else
                % Lower Right -- U
                xrow = row - numbNumCoeffs;
                xcol = col - numbNumCoeffs;
                Q    = 4;

                if (xor(mod(xrow,2) == 0, mod(xcol,2) == 0) == 1) % If both are odd or both are even
                    M(row,col) = 0;
                    power = 0;
                    sign = 0;
%                    fprintf(formatPrint,row,col,'0');
                else
                    power = xrow + xcol;

                    if     (mod(xcol,4) == 1)
                        sign = 1;
                    elseif (mod(xcol,4) == 2) 
                        sign = 1;
                    else
                        sign = -1;
                    end

                    M(row,col) = sign * U(w,real(cData),imag(cData),power);
%                    fprintf(formatPrint,row,col,'U');
                end % else -- check for zeros
%%                fprintf(formatPrint,row,col,xrow,xcol,Q,power,sign);
            end % else -- searching through quadrants
        end % for row = 1:Mrows
    end % for col = 1:Mcols

    % Populate C
    for (row = 1:length(C))
        if (row <= numbNumCoeffs)
            power = row - 1;
            if (mod(power,2) == 0)
                C(row) = S(w,real(cData),power);
            else
                C(row) = T(w,imag(cData),power);
            end
        else
            xrow = row - numbNumCoeffs + 1;
            if (mod(xrow,2) == 0)
                C(row) = 0;
            else
                C(row) = U(w,real(cData),imag(cData),power);
            end
        end
    end
    
    % Solve For N
    N = inv(M) * C;
    M
    N
    C

    numCoeffs = N(1:numbNumCoeffs);
    denCoeffs = [1;N(numbNumCoeffs+1:length(N))]; % Need to add 1 as coeff_0
 
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

