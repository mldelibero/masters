function plotFit(pType, Data1, Data2, x)
    titleSize  = 25;
    legendSize = 20;
    axisSize   = 20;
% Plots different types for the regression analysis
    if pType == plotType.cVectorsDiff
        figure;
        rows = 2;
        cols = 2;

        subplot(rows,cols,1);
        loglog(x,abs(Data1)); hold on;
        loglog(x,abs(Data2));
        title('Magnitude','FontSize',titleSize);
        h_leg = legend('Orig','FitData','Location','north');
        set(h_leg,'FontSize',legendSize);
        xlabel('\omega'  ,'FontSize',axisSize);
        ylabel('Mag (dB)','FontSize',axisSize);

        subplot(rows,cols,2);
        semilogx(x,rad2deg(phase(Data1))); hold on;
        semilogx(x,rad2deg(phase(Data2)));
        title('Phase','FontSize',titleSize);
        h_leg = legend('Orig','FitData','Location','north');
        set(h_leg,'FontSize',legendSize);
        xlabel('\omega','FontSize',axisSize);
        ylabel('\phi (deg)','FontSize',axisSize);

        subplot(rows,cols,3);
        semilogx(x,abs(Data2)-abs(Data1));
        title('Magnitude Error','FontSize',titleSize);
        xlabel('\omega','FontSize',axisSize);
        ylabel('\Delta Mag (\Omega)','FontSize',axisSize);

        subplot(rows,cols,4);
        semilogx(x,rad2deg(phase(Data2))-rad2deg(phase(Data1)));
        title('Phase Error','FontSize',titleSize);
        xlabel('\omega','FontSize',axisSize);
        ylabel('\Delta \phi (deg)','FontSize',axisSize);

    elseif pType == plotType.oneError
        figure;
        rows = 1;
        cols = 1;
        subplot(rows,cols,1);
        plot(Data1);
        title('Norm Magnitude Err^2 + Norm Phase Err^2','FontSize',titleSize);
        xlabel('n','FontSize',axisSize);
        ylabel('Error','FontSize',axisSize);

    elseif pType == plotType.twoErrors
        figure;
        rows = 1;
        cols = 2;

        subplot(rows,cols,1);
        semilogy(Data1);
        title('Magnitude Err^2','FontSize',titleSize);
        xlabel('n','FontSize',axisSize);
        ylabel('Error^2 ((\Delta \Omega)^2)','FontSize',axisSize);

        subplot(rows,cols,2);
        semilogy(Data2);
        title('Phase Error^2','FontSize',titleSize);
        xlabel('n','FontSize',axisSize);
        ylabel('Error^2 ((\Delta \Phi)^2)','FontSize',axisSize);
    end % if pType == plotType.cVectorsDiff
end % function plotFit()

