function fh  = GSEA_Dot_Plot(DATA,nGroups,fWidth,fHight,LegendSizeVal)

FontSize = 7;
minSize = 20;
maxSize = 150;
LineWidth = 0.5;
GridLines = 'on';
% LegendSizeVal = [10 20 30];
% LegendSizeVal = [10 25 50];
% fWidth = 7;
% fHight = 2.5;

if isempty(nGroups)
    nGroups = length(DATA.PATHWAYS.q);
end

xVal = -log10(DATA.PATHWAYS.q(1:nGroups));
xValTxt = '-log_1_0(q-value)';
SizeVal = DATA.PATHWAYS.numGenesInOveralap(1:nGroups);
ColorVal = DATA.PATHWAYS.Ratio(1:nGroups);
ColorlabelTxt = 'GeneRatio';
YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
%YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
YtickLabelTxt = strrep(YtickLabelTxt,'_TARGET_GENES','');
%YtickLabelTxt=strcat('\it',YtickLabelTxt);


Cmap = colorcet('L08');
CLim = [0 ceil((max(ColorVal)+0.001)*100)/100];
Cmap = flipud(Cmap);
Cmap = colormap('winter');
%colormap(cmap);

SizeValPlot = rescale(SizeVal,minSize,maxSize);
LegendSizeValPlot = rescale(LegendSizeVal,minSize,maxSize);

minVal=min(xVal);
maxVal=max(xVal);
rangeVal = maxVal - minVal;
nudgeVal = rangeVal/20;


fh=figure('Name','GSEA Plot','Color','w','Tag','GSEA Plot figure','Units','inches');
fh.Position(3:4) = [fWidth fHight];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,'Units','inches','Clipping','off');
ah.LineWidth = LineWidth;
ah.XGrid = GridLines;
ah.YGrid = GridLines;
ah.YDir = 'Reverse';
ah.YTick = 1:nGroups;
ah.YLim = [0.5 nGroups+0.5];
ah.YTickLabel = YtickLabelTxt;
ah.Colormap = Cmap;
ah.CLim = CLim;
ah.XLim =[minVal-nudgeVal maxVal+nudgeVal];
ah.OuterPosition(3:4) = [fWidth-(fWidth/10) fHight];
ah.TickLength=[ 0.05/nGroups    0.0];
xlabel(ah,xValTxt);

sh = scatter(xVal,1:nGroups,SizeValPlot,ColorVal,'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);

YPos = 1:1.5:1.5*length(LegendSizeValPlot);

shl = scatter(ah.XLim(2)+nudgeVal*1.3,YPos+1,LegendSizeValPlot,[0 0 0]);
text(ah.XLim(2)+2.2*nudgeVal,1,'Count','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)
for i = 1:length(LegendSizeValPlot)
    text(ah.XLim(2)+2.2*nudgeVal,YPos(i)+1,sprintf('%u',LegendSizeVal(i)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
end

ch = colorbar(ah);
ch.Label.String=ColorlabelTxt;
ch.Units="inches";
ch.FontSize=FontSize;
ch.Position=[ah.OuterPosition(1) + ah.OuterPosition(3)+0.12, ah.Position(2) 0.1, fHight/2.5];

