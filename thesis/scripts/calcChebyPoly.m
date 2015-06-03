function [T_k] = calcChebyPoly(k,w)
    % This function calculates a chebyshev polynomial result from a frequency point and an order

    % Normalize a frequency vector to the range of [-1,1]
    if (k < 1)
        T_k = 1;
    elseif (k == 1)
        T_k = w;
    else
        T_1 = 1;
        T_2 = 1;
        T_k = w;

        for (ind = 2:k)
            T_2 = T_1;
            T_1 = T_k;
            T_k = 2*w*T_1 - T_2;
        end
    end % else
end % function [T_k] = calcChebyPoly(k,w)

