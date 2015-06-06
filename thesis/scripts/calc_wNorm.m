function [w_norm] = calc_wNorm(w)
    % This function calculates the Chebyshev polynomials from a set of w

    % Normalize a frequency vector to the range of [-1,1]
    wmin = min(w);
    wmax = max(w);
    for index = 1:length(w)
        %w_norm = 2*(w-w_min) / Cw_max - w_min) - 1
        w_norm(index) = 2*((w(index) - wmin) / (wmax - wmin))-1;
    end
end % function [T,w_norm] = calcChebyPoly(w)
