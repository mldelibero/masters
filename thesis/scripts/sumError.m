function E = sumError(cData,G)
%   Returns the squared error for mag and phase over all of the iterations.
    E = zeros(size(G,1),2);
    for iter = 1:size(G,1)
        for index = 1:size(G,2)
            E(iter,1) = E(iter,1) + (abs  (cData(index))-abs  (G(iter,index)))^2;
            E(iter,2) = E(iter,2) + (phase(cData(index))-phase(G(iter,index)))^2;
        end % for index = 1:size(G,2)
        E(iter,2) = rad2deg(E(iter,2));
    end % for iter = 1:size(G,1)
end % function E = sumError(cData,G)

