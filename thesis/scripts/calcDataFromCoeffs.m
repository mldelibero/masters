function [G,Num,Den] = calcDataFromCoeffs(numCoeff,denCoeff,w)
    % Calculate output from coefficients 
    % G(s) = (a0 + a1*s + a2*s^2 + ... + an*s^n) /
    %        (b0 + b1*s + b2*s^2 + ... + bm*s^m)

    Num = numCoeff(1) * ones(1, length(w));
    Den = denCoeff(1) * ones(1, length(w));

    for h = 2:length(numCoeff);
        Num = Num + numCoeff(h) .* (1i * w).^(h-1);
    end
    for h = 2:length(denCoeff);
        Den = Den + denCoeff(h) .* (1i * w).^(h-1);
    end

    G = Num ./ Den;
end % function G = calcDataFromCoeffs(numCoeff,denCoeff,w)

