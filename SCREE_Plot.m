function fh = SCREE_Plot(PCAModel,Type,FigSize,TitleText)

Type = 'eig';
FontSize = 8;
LineWidth = 1;
MarkerSize = 6;

fh = figure('Name','Bar Plot','Color','w','Tag','Age Scatter Plot',...
    'Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',FontSize,'FontName','Helvetica');


numValues = PCAModel.NumComp;

switch lower(Type)
    case {'eig','eigenvalue'}
        yVal = PCAModel.Eig;
        YLblTxt = "Eigenvalue";
    case {'varcum','vartot','ssxcom','ssxtot','explvarcum'}
        yVal = [0; PCAModel.ExplVarCum];
        YLblTxt = "Cumulative explained variance"
    case {'ssx','explvar'}
        Val = PCAModel.ExplVar;
        YLblTxt = "Cumulative explained variance"


end
    
    plot(ah,yVal,'o','Linewidth',LineWidth,'MarkerSize',MarkerSize,'Color','k')
    hold on
    plot(ah,yVal,'-','Linewidth',LineWidth,'Color','r')
    yline(ah,1,'color','b','Linewidth',LineWidth,'LineStyle','-')
    ah.LineWidth=1;
    ah.XLim = [0.5 length(yVal) + 0.5];
    set(ah,'FontSize',FontSize);
    xlabel('PCA component','FontSize',FontSize,'Interpreter','none')
    ylabel(YLblTxt,'FontSize',FontSize,'Interpreter','none')

   
    title(TitleText,'FontSize',FontSize,'Interpreter','none')
    set(gcf, 'Color', 'w');