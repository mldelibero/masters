function N = regression(filename,modelOrder)
%% Prep
[w, cData] = getData(filename);
numSamples = size(w,1);
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
                
                M(row,col) = sign * lambda(w,power);
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
                M(row,col) = sign * T(w,imag(cData),power);          
            else
                M(row,col) = sign * S(w,real(cData),power);
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
                M(row,col) = sign * S(w,real(cData),power);
            else
                M(row,col) = sign * T(w,imag(cData),power);
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
                
                M(row,col) = sign * U(w,real(cData),imag(cData),power);
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
N = inv(M) * C;

%% Calc Solution
Num = N(1) * ones(size(w,1),1);
for h = 1:modelOrder
    Num = Num + N(h+1) .* (1i * w).^h;
end
Den = ones(size(w,1),1);
for h = 1:modelOrder
    Den = Den + N(h+modelOrder+1) .* (1i * w).^h;
end
G = Num ./ Den;

%% Calculate Residuls
mag_res   = abs(cData) - abs(G);
phase_res = rad2deg(phase(cData)-phase(G));
%% Plot
figure
subplot(2,2,1); 
semilogx(w,abs(G)); hold on; 
semilogx(w,abs(cData));
title('Magnitude');
legend('Model','Emperical');

subplot(2,2,2); 
semilogx(w,rad2deg(phase(G))); hold on; 
semilogx(w,rad2deg(phase(cData)));
title('Phase');

subplot(2,2,3);
semilogx(w,mag_res,'o'); hold on;
title('Residuals of Magnitude');

subplot(2,2,4); 
semilogx(w,phase_res,'x'); hold on; 
title('Residuals of Phase');
%% Sub Functions
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
end
