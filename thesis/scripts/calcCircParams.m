function [params] = calcCircParams(numCoeffs, denCoeffs, modelType)
    % Calculate the circuit parameters for a full model

    params = ones(size(numCoeffs,1), size(numCoeffs,2) + size(denCoeffs,2));

    if modelType == modelTypes.FULL_MODEL
        for index = 1:size(numCoeffs, 1)
            a0 = numCoeffs(index,1);
            a1 = numCoeffs(index,2);
            a2 = numCoeffs(index,3);
            a3 = numCoeffs(index,4);
            b1 = denCoeffs(index,2);
            b2 = denCoeffs(index,3);

            LE(index) = a3 / b2;
            RE(index) = (a2-a3*b1/b2) / b2;
            RL(index) = a0 - RE(index);
            C(index)  = b2*(a0-RE(index)) / (RL(index)*(a1-a3/b2-b1*RE(index)));
            CD(index) = (b1 - b2 / (C(index)*RL(index)) - C(index)*RL(index)) / RL(index);
            RD(index) = b2/(C(index)*CD(index)*RL(index));
        end % for index = 1:size(numCoeffs,2)

        params = [LE', RE', RL', C', CD', RD'];
    end % if modelType == modelTypes.NO_MODEL
end

