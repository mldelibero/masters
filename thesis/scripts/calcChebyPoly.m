function [T,w_norm] = calcChebyPoly(w)
    % This function calculates the Chebyshev polynomials from a set of w

    % Normalize a frequency vector to the range of [-1,1]
    wmin = w(1);
    wmax = w(length(w));
    for index = 1:length(w)
        w_norm(index) = 2*((w(index) - wmin) / (wmax - wmin))-1;
    end
    T(1) = 1;
    T(2) = w_norm(2);
    for index = 3:length(w_norm)
        T(index) = 2 * w_norm(index -1) * T(index - 1) - T(index - 2);
    end
end
