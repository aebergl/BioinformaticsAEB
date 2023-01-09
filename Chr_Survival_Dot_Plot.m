function fh  = Chr_Survival_Dot_Plot(DATA,X_Variable,Y_Samples,SizeType,SizeCutOffs,LegendSizeVal,fWidth,fHight)

FontSize = 8;
minSize = 20;
maxSize = 150;
LineWidth = 0.5;
GridLines = 'on';

% LegendSizeVal = [10 20 30];
% LegendSizeVal = [10 25 50];
% fWidth = 7;
% fHight = 2.5;


if isempty(Y_Samples)
    SampleIndx = ones(DATA.nRow,1);
else 
    SampleIndx = ismember(DATA.RowId,Y_Samples);
end
nGroups = sum (SampleIndx);

indx_X_Val = strcmpi(X_Variable,DATA.ColId);
if any(indx_X_Val)
    xVal = log2(DATA.X(SampleIndx,indx_X_Val));
    xValTxt = 'log_2(HR)';
else
    error('%s not found',X_Variable)
end

% Select Size value
indx_SizeVal = strcmpi(SizeType,DATA.ColId);
if any(indx_SizeVal)
    SizeValOrig = DATA.X(SampleIndx,indx_SizeVal);
    SizeVal = -log10(SizeValOrig);
    %LegendPValsPlot = -log10(SizeCutOffs);
else
    error('%s not found',SizeType)
end


LegendSizeValMat = sum(repmat(SizeCutOffs,nGroups,1) >= repmat(SizeValOrig,1,length(SizeCutOffs)),2);

SizeValPlot = LegendSizeVal(LegendSizeValMat);


%SizeVal = DATA.PATHWAYS.numGenesInOveralap(1:nGroups);
% ColorVal = DATA.PATHWAYS.Ratio(1:nGroups);
% ColorlabelTxt = 'GeneRatio';
% YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
% YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
% YtickLabelTxt = strrep(YtickLabelTxt,'_TARGET_GENES','');
% YtickLabelTxt=strcat('\it',YtickLabelTxt);


% Cmap = colorcet('L08');
% %CLim = [0 ceil((max(ColorVal)+0.001)*100)/100];
% Cmap = flipud(Cmap);
% Cmap = colormap('winter');
%colormap(cmap);

% SizeValPlot = rescale(SizeVal,minSize,maxSize);
% LegendPValsPlot = rescale(LegendPValsPlot,minSize,maxSize);

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
ah.YTickLabel = Y_Samples;
% ah.Colormap = Cmap;
% ah.CLim = CLim;
ah.XLim =[minVal-nudgeVal maxVal+nudgeVal];
ah.OuterPosition(3:4) = [fWidth-(fWidth/10) fHight];
ah.TickLength=[ 0.05/nGroups    0.0];
xlabel(ah,xValTxt);

sh = scatter(xVal,1:nGroups,SizeValPlot,[1 0 0],'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);

YPos = 1:1.5:1.5*length(LegendSizeVal);

shl = scatter(ah.XLim(2)+nudgeVal,YPos+1,LegendSizeVal,[0 0 0]);
text(ah.XLim(2)+1.4*nudgeVal,1,'p-value','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)
SizeCutOffs
for i = 1:length(LegendSizeVal)
    if i==1;
        txt_str = 'N.S.'
    else
        n_val = floor(-log10(SizeCutOffs(i))) + 1;
        txt_str = sprintf('%.*f',n_val,SizeCutOffs(i));
    end
    text(ah.XLim(2)+1.7*nudgeVal,YPos(i)+1,txt_str,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
end

% ch = colorbar(ah);
% ch.Label.String=ColorlabelTxt;
% ch.Units="inches";
% ch.FontSize=FontSize;
% ch.Position=[ah.OuterPosition(1) + ah.OuterPosition(3)+0.02, ah.Position(2) 0.2, fHight/2.5];

