function fh  = Chr_Survival_Dot_Plot(DATA,X_Variable,Y_Samples,SizeType,SizeCutOffs,LegendSizeVal,ColorType,fWidth,fHight)

FontSize = 6;
LineWidth = 0.5;
GridLines = 'on';

% Get samples to use
if isempty(Y_Samples)
    SampleIndx = ones(DATA.nRow,1);
else
    SampleIndx = ismember(DATA.RowId,Y_Samples);
end
nGroups = sum (SampleIndx);

Y_Samples_txt=strrep(Y_Samples,'Hyper ','');


% Select X-value
indx_X_Val = strcmpi(X_Variable,DATA.ColId);
if any(indx_X_Val)
    % xVal = log2(DATA.X(SampleIndx,indx_X_Val));
    % xValTxt = 'log_2(HR)';

    xVal = DATA.X(SampleIndx,indx_X_Val);
    xValTxt = 'log rank HR';
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

% Assign size for each sample
LegendSizeValMat = sum(repmat(SizeCutOffs,nGroups,1) >= repmat(SizeValOrig,1,length(SizeCutOffs)),2);
SizeValPlot = LegendSizeVal(LegendSizeValMat);


% Get colors to use
if isempty(ColorType)
    ColorVal = [0.737254901960784   0.235294117647059   0.160784313725490];
elseif ismatrix(ColorType)
    ColorVal = ColorType;
else
    indx_ColorVal = strcmpi(SizeType,DATA.ColId);
    if any(indx_ColorVal)
        ColorValOrig = DATA.X(SampleIndx,indx_ColorVal);
        ColorVal = -log10(ColorValOrig);
        ColorlabelTxt = '-log_1_0(p-value)';
    else
        error('%s not found',SizeType)
    end

end

% Cmap = colorcet('L08');
CLim = [0 ceil((max(ColorVal)+0.001)*100)/100];
% Cmap = flipud(Cmap);



minVal=min(xVal);
maxVal=max(xVal);
rangeVal = maxVal - minVal;
nudgeVal = rangeVal/10;


fh=figure('Name','GSEA Plot','Color','w','Tag','GSEA Plot figure','Units','inches');
fh.Position(3:4) = [fWidth fHight];
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,'Units','inches','Clipping','off','PositionConstraint','innerposition');
Cmap = colormap('winter');
ah.LineWidth = LineWidth;
ah.XGrid = GridLines;
ah.YGrid = GridLines;
ah.YDir = 'Reverse';
ah.YTick = 1:nGroups;
ah.YLim = [0.5 nGroups+0.5];

%ah.YTickLabel = Y_Samples_txt;
ah.Colormap = Cmap;
ah.CLim = CLim;
ah.XLim =[minVal-nudgeVal maxVal+nudgeVal];
ah.OuterPosition = [0 0 fWidth-0.5 fHight];
ah.TickLength=[ 0.05/nGroups    0.0];
xlabel(ah,xValTxt);

sh = scatter(ah,xVal,1:nGroups,SizeValPlot,ColorVal,'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);

YPos = 2:1.99:1.99*length(LegendSizeVal)+1;

shl = scatter(ah,ah.XLim(2)+nudgeVal,YPos+1,LegendSizeVal,[0 0 0]);
text(ah.XLim(2)+1.4*nudgeVal,1,'p-value','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)

for i = 1:length(LegendSizeVal)
    if i==1;
        txt_str = 'N.S.';
    else
        txt_str = num2str(SizeCutOffs(i));
    end
    text(ah.XLim(2)+1.7*nudgeVal,YPos(i)+1,txt_str,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
end

if ~isempty(ColorType) && ~ismatrix(ColorType)
    ch = colorbar(ah);
    ch.Label.String=ColorlabelTxt;
    ch.Units="inches";
    ch.FontSize=FontSize;
    ch.Position=[fWidth-0.7, ah.Position(2) 0.2, fHight/2.5];
end
