function [Den] = getInitGuess(w,modelType)
    if modelType == modelTypes.FULL_MODEL
        C  = 1;
        RE = 1;
        LE = 1;
        RL = 1;
        CD = 1;
        RD = 1;
        
        b1 = CD*RD+C*RL+CD*RL;
        b2 = C*CD*RD*RL;
        [~,~,Den] = calcDataFromCoeffs([1;1;1],[1,b1,b2],w);

    else % Set denominator to all ones
        Den = ones(length(w));
    end % if modelType == modelTypes.FULL_MODEL
end % function [Den] = getInitGuess(w)
