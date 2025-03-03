function fh = GSEA_DotPlot_MultiDataSets(DATA,XAxisId,XGroups,YAxisId,YGroups,SizeVar,ColorVar,varargin)
% USAGE:
%   fh  = GSEA_DotPlot_MultiDataSets(DATA,nGroups,fWidth,fHight,LegendSizeVal,YTickText)
%   creates a dot plot for GSEA results from multiple datasets
%
% INPUTS:
% * DATA:       DATA structure with GSEA result datasets
% * XAxisId:    Id to be used for groups along the X-axis from RowAnnotationFields.
% * XGroups:    List of Groups to be used and order for the X-axis, [] uses all
% * YAxisId:    Id to be used for groups along the Y-axis from RowAnnotationFields.
% * YGroups:    List of Groups to be used and order for the Y-axis, [] uses all
% * SizeVar:    Name of variable to be used for size from ColId, unless 'VarId' is defined 
% * ColorVar:   Name of variable to be used for color from ColId, unless 'VarId' is defined
%
% OUTPUTS:
%   fh: Figure handle to dot plot figure
%
%   options ---------------------------------------
%
%   'FontSize'      FontSize for all text in figure [7]
%   'FigSize'       Vector with figure width and hight in inches [4 7.5]
%   'LegendSizeVal' Two vectors defining marker size and cut-off values:
%                   Size cut-off for the different sizes [1 0.05 0.01 0.001 0.0001]
%                   Marker size for the different cut-off [5 20 40 60 80]
% 
%   'VarId'         Id for ColAnnotationFields to be used for SizeVar and ColorVar
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


VariableIdentifier = false;
FigSize = [4 7.5];
FontSize = 8;
LineWidth = 0.5;
GridLines = 'on';
RightMargin = 0.7;
SizeCutOffs =[1 0.05 0.01 0.001 0.0001];
LegendSizeVal = [5 20 40 60 80];

CMap = colorcet('D01');
CLim = [-3 3];

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'FigSize')
        i = i + 1;
        FigSize = varargin{i};
    elseif strcmpi(varargin{i},'FontSize')
        i = i + 1;
        FontSize = varargin{i};
    elseif strcmpi(varargin{i},'LegendSizeVal')
        i = i + 1;
        SizeCutOffs = varargin{i};
        i = i + 1;
        LegendSizeVal = varargin{i};
    end
end


% Get groups for X-axis
indx = strcmpi(XAxisId,DATA.RowAnnotationFields);
if ~any(indx)
    error('Error. \n%s not found in DATA.RowAnnotationFields',XAxisId);
elseif sum(indx) > 1
    error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',XAxisId);
else
    XAxisData = DATA.RowAnnotation(:,indx);
end
if isempty(XGroups)
    XAxisData = categorical(XAxisData);
else
    XAxisData = categorical(XAxisData,XGroups);
end


%Get Gene sets
indx = strcmpi(YAxisId,DATA.RowAnnotationFields);
if ~any(indx)
    error('Error. \n%s not found in DATA.RowAnnotationFields',YAxisId);
elseif sum(indx) > 1
    error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',YAxisId);
else
    YAxisData = DATA.RowAnnotation(:,indx);
end
if isempty(YGroups)
    YAxisData = categorical(YAxisData);
else
    YAxisData = categorical(YAxisData,YGroups);
end

% Selection of Size and Colorvariable
if VariableIdentifier
    indx = strcmpi(VariableIdentifier,DATA.ColAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA.ColAnnotationFields',VariableIdentifier);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA.ColAnnotationFields',VariableIdentifier);
    else
        DATA_ID = DATA.ColAnnotation(:,indx);
    end
else
    DATA_ID = DATA.ColId;
end
indx = ismember(DATA_ID,SizeVar);
SizeVal = DATA.X(:,indx);

indx = ismember(DATA_ID,ColorVar);
ColorVal = DATA.X(:,indx);

% Select sample with values
indx_rem = isundefined(XAxisData) | isundefined(YAxisData) | ~isfinite(SizeVal) | ~isfinite(ColorVal);

XAxisData(indx_rem) = [];
YAxisData(indx_rem) = [];
SizeVal(indx_rem) = [];
ColorVal(indx_rem) = [];
nVal = length(ColorVal);
nXval = length(unique(XAxisData));
nYval = length(unique(YAxisData));

% Make numeric for plotting
XAxisDataPlot = double(XAxisData);
XAxisDatalabels = categories(XAxisData);


% Scale according to Legend Size
% SizeValPlot = rescale(SizeVal,minSize,maxSize);
% LegendSizeValPlot = rescale(LegendSizeVal,minSize,maxSize);
% Assign size for each sample
LegendSizeValMat = sum(repmat(SizeCutOffs,nVal,1) >= repmat(SizeVal,1,length(SizeCutOffs)),2);
SizeValPlot = LegendSizeVal(LegendSizeValMat);


fh=figure('Name','GSEA Plot','Color','w','Tag','GSEA Plot','Units','inches');
fh.Position(3:4) = FigSize;
ah = axes(fh,'NextPlot','add','tag','Volcano Plot','box','on','Layer','top','FontSize',FontSize,'Units','inches',...
    'PositionConstraint','outerposition','Clipping','off','YDir','reverse');

ah.LineWidth = LineWidth;
ah.XGrid = GridLines;
ah.YGrid = GridLines;
ah.YDir = 'Reverse';
ah.YTick = 1:nYval;
ah.YLim = [0.5 nYval+0.5];
ah.XTick = 1:nXval;
ah.XLim = [0.5 nXval+0.5];
ah.XAxis.TickLabelRotation=-45;
ah.Colormap = CMap;
ah.CLim = CLim;
ah.OuterPosition(3:4) = [FigSize(1)-RightMargin FigSize(2)];
ah.TickLength=[ 0.05/nVal    0.0];

% switch lower(X_value_name)
%     case 'es'
%         XlabelTxt = "Enrichment score (ES)";
%     case 'nes'
%         XlabelTxt = "Normalized enrichment score (NES)";
% end
% xlabel(ah,XlabelTxt);

sh = scatter(XAxisDataPlot,YAxisData,SizeValPlot,ColorVal,'MarkerFaceColor','flat','MarkerEdgeColor',[0.1 0.1 0.1]);
%ah.XLimMode = 'manual';
% switch lower(GeneSetSource)
%     case 'hallmark'
%                 YtickLabelTxt = string(ah.YAxis.TickLabels);
%                 YtickLabelTxt = strrep(YtickLabelTxt,'HALLMARK_','');
%                 YtickLabelTxt = strrep(YtickLabelTxt,'_',' ');
%                 YtickLabelTxt = lower(YtickLabelTxt);
%                 YtickLabelTxt = upper(extractBefore(YtickLabelTxt,2)) + extractAfter(YtickLabelTxt,1);
% end

% ah.YAxis.TickLabels = YtickLabelTxt;


ah.XAxis.TickLabels = XAxisDatalabels;
% Add color bar
ch = colorbar(ah,'Units','inches','FontSize',FontSize,...
    'Position',[ah.Position(1) + ah.Position(3)+0.1, ah.Position(2) 0.1, FigSize(2)/5]);
ch.Label.String=ColorVar;
ch.FontSize=FontSize;

% Add size legend
StepVal = 1;
nudgeVal = 0.1;
YPos_Size = 1:StepVal:StepVal*length(LegendSizeVal);

X_Start = double(ah.XLim(2)) + 0.5;
shl = scatter(X_Start,YPos_Size,LegendSizeVal,[0 0 0]);
text(X_Start+1.3*nudgeVal,YPos_Size(1)-StepVal,SizeVar,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FontSize)
for i = 1:length(LegendSizeVal)
    if i==1
        txt_str = 'N.S.';
    elseif i==length(LegendSizeVal)
        txt_str = strcat("<",num2str(SizeCutOffs(i)));
    else
        txt_str = num2str(SizeCutOffs(i));
    end
    text(X_Start+2*nudgeVal,YPos_Size(i),txt_str,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',FontSize)
end