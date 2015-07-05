% plotFit.m
% Outputs different plots for the regression analysis.
function [h_fig, h_leg, h_title, ax, p] = plotFit(pType, Data1, Data2, x)
%% Common Prep
titleSize     = 25;
legendSize    = 20;
axisTitleSize = 20;
axisSize      = 15;

h_fig = figure;

%% Plot specific request
if pType == plotType.DIFF_PLOT;
    ax = subplot(1,1,1);
    p(1) = plot(x,Data1,'o'); hold on;
    p(2) = plot(x,Data2);
    h_title = title('Basic LSE','FontSize',titleSize);
    h_leg = legend('Orig','Fit');
    set(h_leg,'FontSize',legendSize);
    xlabel(ax, 'x', 'FontSize', axisTitleSize);
    ylabel(ax, 'y', 'FontSize', axisTitleSize);

elseif pType == plotType.cData
    [ax,p(1),p(2)] = plotyy(x,mag2db(abs(Data1)),x, ...
    rad2deg(phase(Data1)),'semilogx','semilogx');

    h_title = title('GRM31MR71H105KA88','FontSize',titleSize);
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
    colors = [0.000,0.447,0.741 ;
    0.850,0.325,0.098 ;
    0.929,0.694,0.125 ;
    0.494,0.184,0.556 ;
    0.466,0.674,0.188 ;
    0.301,0.745,0.933];

    [ax(1,:),p(1),p(2)] = plotyy(x,mag2db(abs(Data1(1,:))),x,...
    rad2deg(phase(Data1(1,:))),'semilogx','semilogx'); hold on;
    [ax(2,:),p(3),p(4)] = plotyy(x,mag2db(abs(Data1(2,:))),x,...
    rad2deg(phase(Data1(2,:))),'semilogx','semilogx'); hold on;
    [ax(3,:),p(5),p(6)] = plotyy(x,mag2db(abs(Data1(3,:))),x,...
    rad2deg(phase(Data1(3,:))),'semilogx','semilogx');

    h_title = title('RC Filter -- Increasing C','FontSize',titleSize);
    h_leg = legend(p,'Mag1','Pha1','Mag2','Pha2',...
    'Mag3','Pha3','Location','northeast');
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

elseif pType == plotType.MULTDATA
    ax = subplot(1,1,1);

    for ind = 1:size(Data1,1)
        p(ind) = plot(x,Data1(ind,:)); hold on;
    end

    h_title = title('','FontSize',titleSize);
    h_leg = legend();
    set(h_leg,'FontSize',legendSize);
    xlabel(ax(1,1), 'x'  , 'FontSize', axisTitleSize);
    ylabel(ax(1,1), 'y'  , 'FontSize', axisTitleSize);

elseif pType == plotType.cVectorsDiff
    rows = 2;
    cols = 2;

    ax(1) = subplot(rows,cols,1);
    p(1)  = semilogx(x,mag2db(abs(Data1))); hold on;
    p(2)  = semilogx(x,mag2db(abs(Data2)));
    h_title(1) = title('Magnitude','FontSize',titleSize);
    h_leg = legend('Orig','FitData','Location','north');
    set(h_leg,'FontSize',legendSize);
    xlabel('\omega (rad)'  ,'FontSize',axisTitleSize);
    ylabel('Mag (dB)','FontSize',axisTitleSize);
    min_y1 = min(get(ax(1),'ytick'));
    max_y1 = max(get(ax(1),'ytick'));
    set(ax(1),'ytick',min_y1:20:max_y1);

    ax(2) = subplot(rows,cols,2);
    p(3)  = semilogx(x,rad2deg(phase(Data1))); hold on;
    p(4)  = semilogx(x,rad2deg(phase(Data2)));
    h_title(2) = title('Phase','FontSize',titleSize);
    h_leg = legend('Orig','FitData','Location','north');
    set(h_leg,'FontSize',legendSize);
    xlabel('\omega (rad)','FontSize',axisTitleSize);
    ylabel('\phi (deg)','FontSize',axisTitleSize);
    min_y2 = min(get(ax(2),'ytick'));
    max_y2 = max(get(ax(2),'ytick'));
    set(ax(2),'ytick',min_y2:20:max_y2);

    ax(3) = subplot(rows,cols,3);
    p(5)  = semilogx(x,abs(Data2)-abs(Data1));
    h_title(3) = title('Magnitude Error','FontSize',titleSize);
    xlabel('\omega (rad)','FontSize',axisTitleSize);
    ylabel('\Delta Mag (\Omega)','FontSize',axisTitleSize);

    ax(4) = subplot(rows,cols,4);
    p(6)  = semilogx(x,rad2deg(phase(Data2))-rad2deg(phase(Data1)));
    h_title(4) = title('Phase Error','FontSize',titleSize);
    xlabel('\omega (rad)','FontSize',axisTitleSize);
    ylabel('\Delta \phi (deg)','FontSize',axisTitleSize);

elseif pType == plotType.oneError
    rows = 1;
    cols = 1;
    ax = subplot(rows,cols,1);
    p  = plot(Data1);
    h_title = title('Norm Magnitude Err^2 + Norm Phase Err^2','FontSize',titleSize);
    h_leg = legend();
    xlabel('n','FontSize',axisTitleSize);
    ylabel('Error','FontSize',axisTitleSize);

elseif pType == plotType.twoErrors
    rows = 1;
    cols = 2;

    ax(1) = subplot(rows,cols,1);
    p (1) = semilogy(Data1);
    h_title(1) = title('Magnitude Err^2','FontSize',titleSize);
    xlabel('n','FontSize',axisTitleSize);
    ylabel('Error^2 ((\Delta \Omega)^2)','FontSize',axisTitleSize);
    h_leg = legend();

    ax(2) = subplot(rows,cols,2);
    p (2) = semilogy(Data2);
    h_title(2) = title('Phase Error^2','FontSize',titleSize);
    xlabel('n','FontSize',axisTitleSize);
    ylabel('Error^2 ((\Delta \Phi)^2)','FontSize',axisTitleSize);

elseif pType == plotType.MULTPLOT
    rows = 1;
    cols = 2;

    ax(1) = subplot(rows,cols,1);
    p (1) = loglog(x,Data1);
    h_title(1) = title('Magnitude Err^2','FontSize',titleSize);
    xlabel('n','FontSize',axisTitleSize);
    ylabel('Error^2 ((\Delta \Omega)^2)','FontSize',axisTitleSize);
    h_leg = legend();

    ax(2) = subplot(rows,cols,2);
    p (2) = loglog(x,Data2);
    h_title(2) = title('Phase Error^2','FontSize',titleSize);
    xlabel('n','FontSize',axisTitleSize);
    ylabel('Error^2 ((\Delta \Phi)^2)','FontSize',axisTitleSize);

elseif pType == plotType.OPAREA
    [ax,p(1),p(2)] = plotyy(x,Data1,x(3:5),Data2,'semilogx','loglog');

    h_title = title('Safe Operating Area','FontSize',titleSize);
    h_leg = legend('Voltage','Minimum Capacitance','Location','northeast');
    set(h_leg,'FontSize',legendSize);
    xlabel(ax(1), 'R (\Omega)', 'FontSize', axisTitleSize);
    ylabel(ax(1), 'Vmax (V)'  , 'FontSize', axisTitleSize);
    ylabel(ax(2), 'Cmin (uF)' , 'FontSize', axisTitleSize);

    set(ax(1),'ytick',0:100:500);
%    set(ax(2),'YLim',[10^-9 10^-3]);
    ax(2).YTickMode='Auto';
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

