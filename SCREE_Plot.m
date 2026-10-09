function fh = SCREE_Plot(PCAModel,Type,FigSize,TitleText)

if isempty(Type)
    Type = 'eig';
end
FontSize = 8;
LineWidth = 1;
MarkerSize = 6;
InsetStart = 3;

fh = figure('Name','Bar Plot','Color','w','Tag','Age Scatter Plot',...
    'Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Scatter Plot','box','on','Layer','top','FontSize',FontSize,'FontName','Helvetica');


numValues = PCAModel.NumComp;

switch lower(Type)
    case {'eig','eigenvalue'}
        yVal = PCAModel.Eig;
        YLblTxt = "Eigenvalue";
    case {'varcum','vartot','ssxcum','ssxtot','explvarcum'}
        yVal = [0; PCAModel.ExplVarCum];
        YLblTxt = "Cumulative explained variance";
    case {'ssx','explvar'}
        yVal = PCAModel.ExplVar;
        YLblTxt = "Cumulative explained variance";


end

plot(ah,yVal,'o','Linewidth',LineWidth,'MarkerSize',MarkerSize,'Color','k')
hold on
plot(ah,yVal,'-','Linewidth',LineWidth,'Color','r')
yline(ah,1,'color','b','Linewidth',LineWidth,'LineStyle','-')
ah.LineWidth=0.75;
ah.XLim = [0.5 length(yVal) + 0.5];
set(ah,'FontSize',FontSize);
xlabel('PCA component','FontSize',FontSize,'Interpreter','none')
ylabel(YLblTxt,'FontSize',FontSize,'Interpreter','none')
title(TitleText,'FontSize',FontSize,'Interpreter','none')
set(gcf, 'Color', 'w');

yVal_inset = yVal(InsetStart:end);
xVal = InsetStart:length(yVal);
ah_inset = axes(fh,'Position', [0.35, 0.35, 0.5, 0.5],'box','on','Layer','top','FontSize',FontSize,'FontName','Helvetica');
plot(ah_inset,xVal,yVal_inset,'o','Linewidth',LineWidth-0.25,'MarkerSize',MarkerSize,'Color','k')
hold on
plot(ah_inset,xVal,yVal_inset,'-','Linewidth',LineWidth-0.25,'Color','r')
yline(ah_inset,1,'color','b','Linewidth',LineWidth-0.25,'LineStyle','-')
ah_inset.LineWidth=0.75-0.25;
ah_inset.XLim = [xVal(1)-0.5 length(yVal_inset) + 0.5];
set(ah_inset,'FontSize',FontSize-1);
xlabel(ah_inset,'PCA component','FontSize',FontSize-1,'Interpreter','none')
ylabel(ah_inset,YLblTxt,'FontSize',FontSize-1,'Interpreter','none')
