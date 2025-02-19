function fh  = GSEA_Dot_Plot(DATA,nGroups,fWidth,fHight,LegendSizeVal,YTickText,varargin)
% USAGE:
%   fh  = GSEA_Dot_Plot(DATA,nGroups,fWidth,fHight,LegendSizeVal,YTickText)
%   Add column annotation from file to DATA
%
% INPUTS:
% * DATA: DATA structure with results
% * nGroups: Number of entries to display [] shows all.
% * fWidth: Figure width in inches
% * fHight: Figure hight in inches
% * LegendSizeVal: vector with ligend sizees [10 20 30]
% * YTickText: Name variable to be used as size from ColId 'HR coxreg DSS'
% * ColorType: Name variable to be used as size from ColId 'p coxreg DSS'
%
% OUTPUTS:
%   fh: Figure handle to Chromosome figure
%
%   options ---------------------------------------
%
%   'FontSize'      FontSize for all text in figure [7]
%   'REGION'        Highligth one or multiple regions cell structure with 
%                   {[Xstart Xstop], [Ystart Xstop], LabelTxt} 
%   'CytoBand'      Displays the Cytoband nased on HG38, give unit to be used, 'mb'
%   'YValCutOff'    Selects points with abs(value) larger than YValCutOff [0]
%   'SizeLegend'    Displays a size legend defined with the following input vectors:
%                   Size cut-off for the different sizes  [1 0.05 0.01 0.001]
%                   Marker size for the different cut-off [5 10 20 30]%
%                   Y position for the different markers sizes = [2.5 3 3.5 4]
%   'MarkSelected'  Highligh selected points
%   'Print'         Do not transpose the input file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FontSize = 7;
minSize = 20;
maxSize = 100;
LineWidth = 0.5;
GridLines = 'on';
RightMargin = 0.5;
% LegendSizeVal = [10; 20 30];
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

switch DATA.Collection{1}
    case 'H'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                YtickLabelTxt = strrep(YtickLabelTxt,'HALLMARK_','');
                YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
                YtickLabelTxt = lower(YtickLabelTxt);
                YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
            case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end
    case 'CP'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                %YtickLabelTxt = strrep(YtickLabelTxt,'WP_','');
                YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
                YtickLabelTxt = lower(YtickLabelTxt);
                YtickLabelTxt = upper(extractBefore(YtickLabelTxt,' ')) + " " + extractAfter(YtickLabelTxt,' ');
            case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end

    case 'CP:KEGG'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                YtickLabelTxt = strrep(YtickLabelTxt,'KEGG_','');
                YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
                YtickLabelTxt = lower(YtickLabelTxt);
                YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
            case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end
    case 'CP:PID'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                YtickLabelTxt = strrep(YtickLabelTxt,'PID_','');
                YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
                YtickLabelTxt = lower(YtickLabelTxt);
                YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
            case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end

    case 'CP:REACTOME'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                YtickLabelTxt = strrep(YtickLabelTxt,'REACTOME_','');
                YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
                YtickLabelTxt = lower(YtickLabelTxt);
                YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
            case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end
    case 'CP:WIKIPATHWAYS'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                YtickLabelTxt = strrep(YtickLabelTxt,'WP_','');
                YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
                YtickLabelTxt = lower(YtickLabelTxt);
                YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
            case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end
    case 'TFT:GTRD'
        switch lower(YTickText)
            case 'name'
                YtickLabelTxt = DATA.PATHWAYS.Name(1:nGroups);
                YtickLabelTxt = strrep(YtickLabelTxt,'_TARGET_GENES','');
                YtickLabelTxt=strcat('\it',YtickLabelTxt);

           case 'description'
                YtickLabelTxt = DATA.PATHWAYS.Description(1:nGroups);
        end
    
end



Cmap = colorcet('L08');
CLim = [0 ceil((max(ColorVal)+0.001)*100)/100];
Cmap = flipud(Cmap);

%colormap(cmap);

SizeValPlot = rescale(SizeVal,minSize,maxSize);
LegendSizeValPlot = rescale(LegendSizeVal,minSize,maxSize);

minVal=min(xVal);
maxVal=max(xVal);
rangeVal = maxVal - minVal;
nudgeVal = rangeVal/20;


fh=figure('Name','GSEA Plot','Color','w','Tag','GSEA Plot figure','Units','inches');
fh.Position(3:4) = [fWidth fHight];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,'Units','inches',...
    'PositionConstraint','outerposition','Clipping','off');
Cmap = colormap('winter');
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
ah.OuterPosition(3:4) = [fWidth-RightMargin fHight];
ah.TickLength=[ 0.05/nGroups    0.0];
xlabel(ah,xValTxt);

sh = scatter(xVal,1:nGroups,SizeValPlot,ColorVal,'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);

YPos = 1:1.5:1.5*length(LegendSizeValPlot);

shl = scatter(ah.XLim(2)+nudgeVal*1.3,YPos+1,LegendSizeValPlot,[0 0 0]);
text(ah.XLim(2)+2.2*nudgeVal,1,'Count','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)
for i = 1:length(LegendSizeValPlot)
    text(ah.XLim(2)+2.2*nudgeVal,YPos(i)+1,sprintf('%u',LegendSizeVal(i)),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
end

ch = colorbar(ah,'Units','inches','FontSize',FontSize,...
    'Position',[ah.Position(1) + ah.Position(3)+0.1, ah.Position(2) 0.1, fHight/2.5]);
ch.Label.String=ColorlabelTxt;
ch.FontSize=FontSize;
%ch.Position=[fWidth-RightMargin-0.4, ah.Position(2) 0.1, fHight/2.5];

