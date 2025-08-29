function [fh, varargout]= VolcanoPlotResults(DATA,X_Variable,X_CutOff,Y_Variable,Y_CutOff,varargin)
AlphaRange = [0.01 0.5];
printResults = false;
FontSize = 8;
EqualXLim = false;
PlotType = 'simple';
Delimiter = '\t';
SizeRange = [1 75];
FigureSize = [1.8 1.8];
TopNValues = 0;
TopPrctile = 0;
RowsToHighligh = [];
Cmap_UpDn = [1 0 0; 0 0 1];
XLim = [];
YLim = [];
TitleTxt = false;
LineWidth = 0.5;
DataTipVariableName = [];
XlimCrop = false;

MarkerLineWidth = 0.1;
MarkerEdgeColor = [0 0 0];
MarkerLineAlphaValue = 0.8;


i=0;
%Check input options
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'AlphaRange')
        i = i + 1;
        AlphaRange = varargin{i};
    elseif strcmpi(varargin{i},'SizeRange')
        i = i + 1;
        SizeRange = varargin{i};
    elseif strcmpi(varargin{i},'FigureSize')
        i = i + 1;
        FigureSize = varargin{i};
    elseif strcmpi(varargin{i},'TopNValues')
        i = i + 1;
        TopNValues = varargin{i};
    elseif strcmpi(varargin{i},'TopPrctile')
        i = i + 1;
        TopPrctile = varargin{i};
    elseif strcmpi(varargin{i},'RowsToHighligh')
        i = i + 1;
        RowsToHighligh = varargin{i};
    elseif strcmpi(varargin{i},'UpDnColors')
        i = i + 1;
        Cmap_UpDn = varargin{i};
    elseif strcmpi(varargin{i},'XLim')
        i = i + 1;
        XLim = varargin{i};
    elseif strcmpi(varargin{i},'YLim')
        i = i + 1;
        YLim = varargin{i};
    elseif strcmpi(varargin{i},'EqualXLim')
        EqualXLim = true;
    elseif strcmpi(varargin{i},'Print')
        printResults = true;    
    elseif strcmpi(varargin{i},'XlimCrop')
        XlimCrop = true;
    elseif strcmpi(varargin{i},'FontSize')
        i = i + 1;
        FontSize = varargin{i};
    elseif strcmpi(varargin{i},'Title')
        i = i + 1;
        TitleTxt = varargin{i};
    elseif strcmpi(varargin{i},'DataTip')
        i = i + 1;
        DataTipVariableName = varargin{i};
    else 
        error('%s is not a valid input option',varargin{i}')
    end
end


% Select X value
indx_X_Val = strcmpi(X_Variable,DATA.ColId);
if any(indx_X_Val)
    x_data = DATA.X(:,indx_X_Val);
else
    error('%s not found',X_Variable)
end

% Not a perfect solution but it works
switch X_Variable
    case 'Delta Average'
        if strcmp('M-value',DATA.Info.DataType)
            XLabel = {'\Delta M-value'};
        elseif strcmp('beta-value',DATA.Info.DataType)
            XLabel = {'\Delta \beta-value'};
        else
        XLabel = {'Log_2 FC'};
        end

    case 'Delta Median'
        XLabel = {'\Delta \beta-value'};
    case 'HR logrank DSS'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR DSS)'};
    case 'HR coxreg DSS'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR DSS)'};
    case 'HR logrank PFS'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR PFS)'};
    case 'HR coxreg PFS'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR PFS)'};
    case 'HR logrank PFI'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR PFI)'};
    case 'HR coxreg PFI'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR PFI)'};
    case 'HR logrank OS'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR OS)'};
    case 'HR coxreg OS'
        x_data =  log2(x_data);
        XLabel = {'log_2(HR OS)'};
    case 'r Spearman'
        XLabel = {'Spearman''s \rho'};
    case 'log2FoldChange'
        XLabel = {'Log_2 FC'};

end

% Select Y value
indx_Y_Val = strcmpi(Y_Variable,DATA.ColId);
if any(indx_Y_Val)
    y_data = DATA.X(:,indx_Y_Val);
else
    error('%s not found',Y_Variable)
end
indx_Zero = y_data == 0;
minVal = min(y_data(~indx_Zero));
y_data(indx_Zero) = minVal;


% Not a perfect solution but it works
switch Y_Variable
    case 'p logrank PFS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p logrank PFI)'};
    case 'q logrank PFS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q logrank PFI)'};
    case 'fdr logrank PFS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr logrank PFI)'};
    case 'p coxreg PFS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p coxreg PFI)'};
    case 'q coxreg PFS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q coxreg PFI)'};
    case 'fdr coxreg PFS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr coxreg PFI)'};
    case 'p logrank PFI'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p logrank PFI)'};
    case 'q logrank PFI'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q logrank PFI)'};
    case 'fdr logrank PFI'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr logrank PFI)'};
    case 'p coxreg PFI'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p coxreg PFI)'};
    case 'q coxreg PFI'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q coxreg PFI)'};
    case 'fdr coxreg PFI'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr coxreg PFI)'};
    case 'p logrank DSS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p logrank DSS)'};
    case 'q logrank DSS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q logrank DSS)'};
    case 'fdr logrank DSS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr logrank DSS)'};
    case 'p coxreg DSS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p coxreg DSS)'};
    case 'q coxreg DSS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q coxreg DSS)'};
    case 'fdr coxreg DSS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr coxreg DSS)'};
    case 'p logrank OS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p logrank OS)'};
    case 'q logrank OS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q logrank OS)'};
    case 'fdr logrank OS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr logrank OS)'};
    case 'p coxreg OS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p coxreg OS)'};
    case 'q coxreg OS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q coxreg OS)'};
    case 'fdr coxreg OS'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr coxreg OS)'};
    case 'p t-test'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p t-test)'};
    case 'q t-test'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q t-test)'};
    case 'fdr t-test'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr t-test)'};
    case 'p MW-test'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p MW-test)'};
    case 'q MW-test'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(q MW-test)'};
    case 'fdr MW-test'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(fdr MW-test)'};
    case 'p Spearman'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p Spearman)'};
    case 'padj'
        y_data =  -log10(y_data);
        YLabel = {'-log_1_0(p-adjusted)'};

end

%Create Figure
fh=figure('Name','Volcano Plot','Color','w','Tag','Volcano Plot figure','Units','inches');
fh.Position(3:4) = FigureSize;
%Create Axis
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','off','Layer','top','FontSize',FontSize,'YAxisLocation','origin','PositionConstraint','outerposition');

ah.LineWidth = LineWidth;
switch PlotType
    case 'simple'
        ah.XGrid= 'off';
        ah.YGrid= 'off';
    otherwise
        ah.XGrid= 'on';
        ah.YGrid= 'on';
end

if TopNValues
    Y_CutOff = min(maxk(y_data,TopNValues+1));
end

if TopPrctile
    Y_CutOff = prctile(y_data,TopPrctile);
end
MaxVal_X_Significant = max(abs(x_data(y_data > Y_CutOff)));

if XlimCrop
    x_data(x_data>MaxVal_X_Significant) = MaxVal_X_Significant;
    x_data(x_data<-MaxVal_X_Significant) = -MaxVal_X_Significant;
end

dist_val = pdist2([x_data y_data],[0 0],"seuclidean");
dist_alpha = rescale(dist_val.^2,AlphaRange(1),AlphaRange(2));
dist_size = rescale(dist_val.^2,SizeRange(1),SizeRange(2));

indx_pos = x_data > 0;
indx_neg = x_data < 0;

x_data_pos = x_data(indx_pos);
y_data_pos = y_data(indx_pos);
dist_pos_size = dist_size(indx_pos);
dist_pos_alpha = dist_alpha(indx_pos);

x_data_neg = x_data(indx_neg);
y_data_neg = y_data(indx_neg);
dist_neg_size = dist_size(indx_neg);
dist_neg_alpha = dist_alpha(indx_neg);

indx_pos_scatter = x_data_pos > X_CutOff & y_data_pos > Y_CutOff;
indx_neg_scatter = x_data_neg < -X_CutOff & y_data_neg > Y_CutOff;

% cMap=colormap('bone');
% cMap=flipud(cMap);
cMap = colormap(colorcet('L02','reverse',true));

% Select Sample Id for data tip
if isempty(DataTipVariableName)
    SampleId = DATA.RowId;
else
    indx_DataTip = strcmpi(DataTipVariableName,DATA.RowAnnotationFields);
    SampleId = DATA.RowAnnotation(:,indx_DataTip);
end

SampleId_pos = SampleId(indx_pos);
SampleId_neg = SampleId(indx_neg);

if any(indx_pos_scatter)
    DensScat(x_data_pos(~indx_pos_scatter),y_data_pos(~indx_pos_scatter),'TargetAxes',ah,'ColorBar',false,'ColorMap',cMap,'mSize',25,'AxisType','auto');
    sh_pos = scatter(ah,x_data_pos(indx_pos_scatter),y_data_pos(indx_pos_scatter),dist_pos_size(indx_pos_scatter),Cmap_UpDn(1,:),'MarkerFaceColor','flat','MarkerEdgeColor',MarkerEdgeColor);
    sh_pos.AlphaDataMapping = 'none';
    sh_pos.AlphaData = dist_pos_alpha(indx_pos_scatter);
    sh_pos.MarkerFaceAlpha = 'flat';
    sh_pos.MarkerEdgeAlpha = MarkerLineAlphaValue;
    sh_pos.LineWidth = MarkerLineWidth;

    row = dataTipTextRow('',SampleId_pos(indx_pos_scatter));
    sh_pos.DataTipTemplate.DataTipRows = row;
    sh_pos.DataTipTemplate.FontSize = FontSize;
    sh_pos.DataTipTemplate.FontAngle = 'normal';
    sh_pos.DataTipTemplate.Interpreter = 'none';
end

if any(indx_neg_scatter)
    DensScat(x_data_neg(~indx_neg_scatter),y_data_neg(~indx_neg_scatter),'TargetAxes',ah,'ColorBar',false,'ColorMap',cMap,'mSize',25,'AxisType','auto');
    sh_neg = scatter(ah,x_data_neg(indx_neg_scatter),y_data_neg(indx_neg_scatter),dist_neg_size(indx_neg_scatter),Cmap_UpDn(2,:),'MarkerFaceColor','flat','MarkerEdgeColor',MarkerEdgeColor);
    sh_neg.AlphaDataMapping = 'none';
    sh_neg.AlphaData = dist_neg_alpha(indx_neg_scatter);
    sh_neg.MarkerFaceAlpha = 'flat';
    sh_neg.MarkerEdgeAlpha = MarkerLineAlphaValue;
    sh_neg.LineWidth = MarkerLineWidth;

    row = dataTipTextRow('',SampleId_neg(indx_neg_scatter));
    sh_neg.DataTipTemplate.DataTipRows = row;
    sh_neg.DataTipTemplate.FontSize = FontSize;
    sh_neg.DataTipTemplate.FontAngle = 'normal';
    sh_neg.DataTipTemplate.Interpreter = 'none';
end

if ~isempty(RowsToHighligh)
    indx = ismember(DATA.RowId,RowsToHighligh);

    scatter(ah,x_data(indx),y_data(indx),3,[0 0 0],'filled');

end




xlabel(ah,XLabel);
ylabel(ah,YLabel);

min_x=min(x_data_neg,[],'all','omitnan');
max_x=max(x_data_pos,[],'all','omitnan');
nudge_x = (max_x-min_x)/20;


if isempty(XLim)
    ah.XLim = [min_x-nudge_x max_x+nudge_x];
else
    ah.XLim = XLim;
end

if EqualXLim
    ah.XLim = [-max(ah.XLim), max(ah.XLim)];
end

min_y=0;
max_y=max(y_data,[],'all','omitnan');
nudge_y = (max_y-min_y)/20;

if isempty(YLim)
    ah.YLim=[0 max_y+nudge_y];
else
    ah.YLim = YLim;
end

ah.YAxis.Label.HorizontalAlignment='center';
ah.YAxis.Label.VerticalAlignment='bottom';

if TitleTxt
    th = title(TitleTxt,'FontWeight','normal','FontSize',FontSize+2);
    th.Position(2) =th.Position(2) + nudge_y*1.5;
    ah.Position(4) = 0.75;
end


switch PlotType
    case 'simple'
        pos_str=sprintf('n=%u',sum(indx_pos_scatter));
        text(ah,ah.XLim(2),ah.YLim(2),pos_str,'Clipping','off','FontSize',FontSize,'HorizontalAlignment','right','VerticalAlignment','top' );

        neg_str=sprintf('n=%u',sum(indx_neg_scatter));
        text(ah,ah.XLim(1),ah.YLim(2),neg_str,'Clipping','off','FontSize',FontSize,'HorizontalAlignment','left','VerticalAlignment','top' );

    otherwise
        for i = 2:length(ah.YAxis.TickValues)
            TickVal = ah.YAxis.TickValues(i);
            nPos = sum(y_data_pos>TickVal);
            pos_str=sprintf('n=%u',nPos);
            text(ah,ah.XLim(2)-nudge_x,TickVal,pos_str,'Clipping','off','HorizontalAlignment','right','FontSize',FontSize )

            nNeg = sum(y_data_neg>TickVal);
            neg_str=sprintf('n=%u',nNeg);
            text(ah,ah.XLim(1)+nudge_x,TickVal,neg_str,'Clipping','off','HorizontalAlignment','left','FontSize',FontSize)
        end
end
ah.Position([1 3]) = [0.15 0.7];

    indx  = find(abs(x_data) > X_CutOff & y_data > Y_CutOff);

if nargout == 2
    varargout{1} = DATA.RowId(indx);

end

if printResults
    format_str_txt = sprintf('%s%%s',Delimiter);
    format_str_val = sprintf('%s%%g',Delimiter);
    format_str_short = sprintf('%%s');
    fprintf('Id');
    fprintf(format_str_txt,DATA.RowAnnotationFields{:});
    fprintf(format_str_txt,DATA.ColId{:});
    fprintf('\n');
    for i=1:length(indx)
        fprintf('%s',DATA.RowId{indx(i)});
        fprintf(format_str_txt,DATA.RowAnnotation{indx(i),:});
        fprintf(format_str_val,DATA.X(indx(i),:));
        fprintf('\n');
    end
end
