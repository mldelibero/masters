function plotFit(pType, Data1, Data2, x)
% Plots different types for the regression analysis
    %% Common Prep
    titleSize     = 25;
    legendSize    = 20;
    axisTitleSize = 20;
    axisSize      = 15;
    
    f_handle = figure;

    %% Plot specific request
    if      pType == plotType.DIFF_PLOT;
        ax = subplot(1,1,1);
        p(1) = plot(x,Data1,'o'); hold on;
        p(2) = plot(x,Data2);
        title('Basic LSE','FontSize',titleSize);
        h_leg = legend('Orig','Fit');
        set(h_leg,'FontSize',legendSize);
        xlabel(ax, 'x', 'FontSize', axisTitleSize);
        ylabel(ax, 'y', 'FontSize', axisTitleSize);

    elseif pType == plotType.cData
        [ax,p(1),p(2)] = plotyy(x,mag2db(abs(Data1)),x,rad2deg(phase(Data1)),'semilogx','semilogx');

        title('GRM31MR71H105KA88','FontSize',titleSize);
        h_leg = legend('Magnitude','Phase','Location','north');
        set(h_leg,'FontSize',legendSize);
        xlabel(ax(1), '\omega (rad)'  , 'FontSize', axisTitleSize);
        ylabel(ax(1), 'Mag (dB)'      , 'FontSize', axisTitleSize);
        ylabel(ax(2), '\phi (deg)'    , 'FontSize', axisTitleSize);
        min_y1 = min(get(ax(1),'ytick'));
        max_y1 = max(get(ax(1),'ytick'));
        min_y2 = min(get(ax(2),'ytick'));
        max_y2 = max(get(ax(2),'ytick'));
        set(ax(1),'ytick',min_y1:20:max_y1);
        set(ax(2),'ytick',min_y2:20:max_y2);

    elseif pType == plotType.MULTCDATA
        colors = [0                ,0.447000000000000,0.741000000000000 ;
                  0.850000000000000,0.325000000000000,0.0980000000000000;
                  0.929000000000000,0.694000000000000,0.125000000000000 ;
                  0.494000000000000,0.184000000000000,0.556000000000000 ;
                  0.466000000000000,0.674000000000000,0.188000000000000 ;
                  0.301000000000000,0.745000000000000,0.933000000000000];
        
        [ax(1,:),p(1),p(2)] = plotyy(x,mag2db(abs(Data1(1,:))),x,rad2deg(phase(Data1(1,:))),'semilogx','semilogx'); hold on;
        [ax(2,:),p(3),p(4)] = plotyy(x,mag2db(abs(Data1(2,:))),x,rad2deg(phase(Data1(2,:))),'semilogx','semilogx'); hold on;
        [ax(3,:),p(5),p(6)] = plotyy(x,mag2db(abs(Data1(3,:))),x,rad2deg(phase(Data1(3,:))),'semilogx','semilogx');

        title('RC Filter -- Increasing C','FontSize',titleSize);
        h_leg = legend(p,'Mag1','Pha1','Mag2','Pha2','Mag3','Pha3','Location','northeast');
        set(h_leg,'FontSize',legendSize);
        xlabel(ax(1,1), '\omega (rad)'  , 'FontSize', axisTitleSize);
        ylabel(ax(1,1), 'Mag (dB)'      , 'FontSize', axisTitleSize);
        ylabel(ax(1,2), '\phi (deg)'    , 'FontSize', axisTitleSize);

        for ind = 1:6
            p(ind).Color = colors(ind,:);
            ax(ind).YColor = [0,0,0]; % Black
        end

        p(2).Marker = 'o';
        p(4).Marker = 'o';
        p(6).Marker = 'o';

    elseif pType == plotType.cVectorsDiff
        rows = 2;
        cols = 2;

        ax(1) = subplot(rows,cols,1);
        p(1)  = semilogx(x,mag2db(abs(Data1))); hold on;
        p(2)  = semilogx(x,mag2db(abs(Data2)));
        title('Magnitude','FontSize',titleSize);
        h_leg = legend('Orig','FitData','Location','north');
        set(h_leg,'FontSize',legendSize);
        xlabel('\omega'  ,'FontSize',axisTitleSize);
        ylabel('Mag (dB)','FontSize',axisTitleSize);
        min_y1 = min(get(ax(1),'ytick'));
        max_y1 = max(get(ax(1),'ytick'));
        set(ax(1),'ytick',min_y1:20:max_y1);

        ax(2) = subplot(rows,cols,2);
        p(3)  = semilogx(x,rad2deg(phase(Data1))); hold on;
        p(4)  = semilogx(x,rad2deg(phase(Data2)));
        title('Phase','FontSize',titleSize);
        h_leg = legend('Orig','FitData','Location','north');
        set(h_leg,'FontSize',legendSize);
        xlabel('\omega','FontSize',axisTitleSize);
        ylabel('\phi (deg)','FontSize',axisTitleSize);
        min_y2 = min(get(ax(2),'ytick'));
        max_y2 = max(get(ax(2),'ytick'));
        set(ax(2),'ytick',min_y2:20:max_y2);

        ax(3) = subplot(rows,cols,3);
        p(5)  = semilogx(x,abs(Data2)-abs(Data1));
        title('Magnitude Error','FontSize',titleSize);
        xlabel('\omega','FontSize',axisTitleSize);
        ylabel('\Delta Mag (\Omega)','FontSize',axisTitleSize);

        ax(4) = subplot(rows,cols,4);
        p(6)  = semilogx(x,rad2deg(phase(Data2))-rad2deg(phase(Data1)));
        title('Phase Error','FontSize',titleSize);
        xlabel('\omega','FontSize',axisTitleSize);
        ylabel('\Delta \phi (deg)','FontSize',axisTitleSize);

    elseif pType == plotType.oneError
        rows = 1;
        cols = 1;
        ax = subplot(rows,cols,1);
        p  = plot(Data1);
        title('Norm Magnitude Err^2 + Norm Phase Err^2','FontSize',titleSize);
        xlabel('n','FontSize',axisTitleSize);
        ylabel('Error','FontSize',axisTitleSize);

    elseif pType == plotType.twoErrors
        rows = 1;
        cols = 2;

        ax(1) = subplot(rows,cols,1);
        p (1) = semilogy(Data1);
        title('Magnitude Err^2','FontSize',titleSize);
        xlabel('n','FontSize',axisTitleSize);
        ylabel('Error^2 ((\Delta \Omega)^2)','FontSize',axisTitleSize);

        ax(2) = subplot(rows,cols,2);
        p (2) = semilogy(Data2);
        title('Phase Error^2','FontSize',titleSize);
        xlabel('n','FontSize',axisTitleSize);
        ylabel('Error^2 ((\Delta \Phi)^2)','FontSize',axisTitleSize);
    end % if pType == plotType.cVectorsDiff

    %% Common plotting options
    for ind = 1:length(ax)
        grid(ax(ind),'on');
        set(ax(ind),'FontSize',axisSize);
    end % for ind = 1:length(ax)

    for ind = 1:length(p)
        p(ind).LineWidth = 2;
    end % for ind = 1:length(p)
end % function plotFit()

    
