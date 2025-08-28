function fh = PlotBoxPlotDATA(DATA,VariableId,GroupVariableName,GroupsToUse,varargin)
% USAGE:
%   fh  = PlotBoxPlotDATA(DATA,VariableId,GroupVariableName,GroupsToUse,varargin)
%   creates a scatter boxplot plot from DATA structure
%
% INPUTS:
% * DATA:       DATA structure w
% * VariableId: Id to be used for groups along the X-axis from RowAnnotationFields.
% * GroupVariableName:    List of Groups to be used and order for the X-axis, [] uses all
% * GroupsToUse:    Id to be used for groups along the Y-axis from RowAnnotationFields.
%
% OUTPUTS:
%   fh: Figure handle to dot plot figure
%
%   options ---------------------------------------
%
%   'FontSize'          FontSize for all text in figure [7]
%   'FigSize'           Vector with figure width and hight in inches [ [2 3.5] ]
%   'MarkerSize'        Marker size for scatter points [30]
%   'MarkerLineWidth'   Marker line width for scatter points [1]
%   'MarkerType'        Marker types for scatter points,{'o','d','^'} [{'o'}]
%   'MarkerFaceColor'   Marker edge color for scatter points ['none']
%   'ColorMap'          Marker edge color for scatter points [GetPalette('Science')]
%   'XJitterWidth'      Width for scatter points [0.6]
%   'BoxLineWidth'      LineWidth for boxes [1]
%   'BoxColor'          Color for box lines [ [0 0 0] ]
%   'BoxWidth'          Width for boxes [0.8]
%   'YlabelTxt'         Text to be used for Y label, if not defined VariableId will be used
%   'TitleTxt'          text for Title
%   'CalcStats'         Calculate stats between groups, Nx2 matrix defines comparisons, [] all
%   'StatType'          Type of comparison, Mann Whitney (MW) or T-test (t-test) ['MW']
%   'PlotStars'         Use starts instead of p-values
%   'StatLineWidth'     Line witth for line btween groups [0.5]
%   'Show_NS'           Show not significant results [false]
%   'TargetAxes'        Axes handle for target axist []
%   'XTickAngle'        Angle for X ticks [-45]
%   'SortData'          Sort group order based on median value
%   'MultipleY'         How to summarize multiple y's Mean' or 'PCA' ['Mean']
%   'DataTipId'         Id to be used for datatip 'RowId' or id from 'RowAnnotationFields' ['RowId']

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FontSize = 10;
FigSize = [2,3.5];

MarkerSize = 30;
MarkerLineWidth = 1;
MarkerType = {'o'};
MarkerFaceColor = 'none';
ColorMap = GetPalette('Science');
XJitterWidth = 0.6;

BoxLineWidth = 1;
BoxColor = [0 0 0];
BoxWidth = 0.8;

YlabelTxt = [];
TitleTxt = [];
VariableIdentifier = false;

CalcStats = false;
CalcGroup = [];
CalcGroupAllUnique = true;
StatType = 'MW';
PlotStars = true;
StatLineWidth = 0.5;
Show_NS=true;

TargetAxes = false;
XTickAngle = -45;
SortData=false;
MultipleY = 'Mean';
DataTipId = 'RowId';


% Check input
if nargin > 4
    ArgsList = {'VariableIdentifier','FigSize','YlabelTxt','TitleText','ColorMap','MarkerTypes',...
        'CalcStats','TargetAxes','MarkerSize','MarkerLineWidth','BoxLineWidth','FontSize',...
        'BoxWidth','XJitterWidth','StatType','PlotStars','Show_NS','XTickAngle','SortData',...
        'MultipleY','DataTipId','BoxColor','MarkerType','MarkerFaceColor'};
    for j=1:2:numel(varargin)

        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)

                case 'variableidentifier'
                    VariableIdentifier = ArgVal;
                case 'figsize'
                    FigSize = ArgVal;
                case 'ylabel'
                    YlabelTxt =ArgVal;
                case 'titletext'
                    TitleTxt = ArgVal;
                case 'colormap'
                    ColorMap = ArgVal;
                case 'markersize'
                    MarkerSize = ArgVal;
                case 'markerlinewidth'
                    MarkerLineWidth = ArgVal;
                case 'markertypes'
                    MarkerType = ArgVal;
                case 'markerfacecolor'
                    MarkerFaceColor = ArgVal;
                case 'calcstats'
                    CalcGroup = ArgVal;
                    CalcStats = true;
                case 'targetaxes'
                    TargetAxes = ArgVal;
                case 'boxlinewidth'
                    BoxLineWidth = ArgVal;
                case 'boxcolor'
                    BoxColor = ArgVal;
                case 'fontsize'
                    FontSize = ArgVal;
                case 'boxwidth'
                    BoxWidth = ArgVal;
                case 'xjitterwidth'
                    XJitterWidth = ArgVal;
                case 'stattype'
                    StatType = ArgVal;
                case 'plotstars'
                    PlotStars = ArgVal;
                case 'show_ns'
                    Show_NS = true;
                case 'xtickangle'
                    XTickAngle = ArgVal;
                case 'sortdata'
                    SortData = ArgVal;
                case 'multipley'
                    MultipleY = ArgVal;
                case 'datatipid'
                    DataTipId = ArgVal;
            end

        end
    end
end

% Select Samples
if isempty(GroupVariableName)
    nGroups = 1;
    SampleIndxMat = true(DATA.nRow,nGroups);
    GroupName = [];
else
    indx = strcmpi(GroupVariableName,DATA.RowAnnotationFields);
    if ~any(indx)
        error('Error. \n%s not found in DATA.RowAnnotationFields',SampleIdentifier);
    elseif sum(indx) > 1
        error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',SampleIdentifier);
    else
        GroupData = DATA.RowAnnotation(:,indx);
    end

    if isempty(GroupsToUse)
        GroupsToUse = unique(GroupData,'stable');
    end
    nGroups = length(GroupsToUse);
    numValuesPerGroup = zeros(nGroups,1);
    SampleIndxMat = false(DATA.nRow,nGroups);
    GroupVariableName = cell(DATA.nRow,1);
    GroupVariableNumber = ones(DATA.nRow,1) * NaN;
    GroupName = cell(nGroups,1);
    GroupNumber = zeros(nGroups,1);
    for i = 1:nGroups
        tmp_name = GroupsToUse{i};
        if ~iscell(tmp_name)
            tmp_name = {tmp_name};
        end
        SampleIndxMatrix = false(DATA.nRow,numel(tmp_name));
        for j=1:length(tmp_name)
            SampleIndxMatrix(:,j) = strcmp(tmp_name{j},GroupData);
        end
        SampleIndxMat(:,i) = any(SampleIndxMatrix,2);
        GroupVariableName(SampleIndxMat(:,i)) = tmp_name(1);
        GroupVariableNumber(SampleIndxMat(:,i)) = i;
        numValuesPerGroup(i) = sum(SampleIndxMat(:,i));
        GroupName(i) = tmp_name(1);
        GroupNumber(i) = i;
    end
end

if isempty(CalcGroup)
    CalcGroup = ones(nGroups*(nGroups-1)/2,2);
    counter = 0;
    CalcGroupAllUnique = false;
    for i=1:nGroups - 1
        for j=i+1:nGroups
            counter = counter + 1;
            CalcGroup(counter,:) = [ i j];

        end
    end

end

SampleIndxToUse = any(SampleIndxMat,2);

% Selection of Y variable
if VariableIdentifier
    indx_VarId = strcmpi(VariableIdentifier,DATA.ColAnnotationFields);
    if ~any(indx_VarId)
        error('Error. \n%s not found in DATA.ColAnnotationFields',VariableIdentifier);
    elseif sum(indx_VarId) > 1
        error('Warning. \nMultiple matches for %s found in DATA.ColAnnotationFields',VariableIdentifier);
    else
        DATA_ID = DATA.ColAnnotation(:,indx_VarId);
    end
else
    DATA_ID = DATA.ColId;
end
indx = ismember(DATA_ID,VariableId);
y_var = DATA.X(:,indx);
y_var(~SampleIndxToUse) = NaN;


if size(y_var,2) > 1
    switch lower(MultipleY)
        case 'mean'
            y_var=mean(y_var,2,"omitnan");
            VariableId = "Average";
        case 'median'
            y_var=median(y_var,2,"omitnan");
            VariableId = "Median";

        case 'pca'
            [PCA_Model] = NIPALS_PCA(y_var,'NumComp',2,'ScaleX',false);
            y_var = PCA_Model.T(:,1);
            VariableId = "PC1";
    end
end

% ColorMap
ColorMap = repmat(ColorMap,ceil(nGroups/size(ColorMap,1)),1);
ColorMap = ColorMap(1:nGroups,:);

% Marker Type
if size(MarkerType,1) < size(MarkerType,2)
    MarkerType = MarkerType';
end
MarkerType = repmat(MarkerType,ceil(nGroups/size(MarkerType,1)),1);
MarkerType = MarkerType(1:nGroups,:);

% sort data
if SortData
    MedianX = zeros(nGroups,1);
    for i=1:nGroups
        indx_sort = SampleIndxMat(:,i);
        MedianX(i) = median(y_var(indx_sort),'omitnan');
    end
    [~,indx_sort] = sort(MedianX,SortData);
    GroupVariableNumberNew = GroupVariableNumber;
    for i=1:nGroups
        GroupVariableNumberNew(GroupVariableNumber == indx_sort(i)) = i;
    end
    GroupVariableNumber= GroupVariableNumberNew;
    SampleIndxMat = SampleIndxMat(:,indx_sort);
    GroupName = GroupName(indx_sort);
end

% Create Figure
if isgraphics(TargetAxes,'axes')
    ah = TargetAxes;
    fh = TargetAxes.Parent;
else
    fh = figure('Name','Box Plot','Color','w','Tag','Box Plot','Units','inches','Colormap',ColorMap);
    fh.Position(3:4) = FigSize;
    fh.Renderer='painters';
    ah = axes(fh,'NextPlot','add','tag','Gene Sample Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
        'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');
    ah.LineWidth = 0.5;
end
ah.Colormap=ColorMap;

% Select Sample Id for data tip text

switch lower(DataTipId)
    case 'rowid'
        SampleId = DATA.RowId;
    otherwise
        indx_SampleId = strcmpi(DataTipId,DATA.RowAnnotationFields);
        if ~any(indx_SampleId)
            error('Error. \n%s not found in DATA.RowAnnotationFields',DataTipId);
        elseif sum(indx_VarId) > 1
            error('Warning. \nMultiple matches for %s found in DATA.RowAnnotationFields',DataTipId);
        else
            SampleId = DATA.RowAnnotation(:,indx_SampleId);
        end
end

% Create scatter points
for i=1:nGroups
    indx = SampleIndxMat(:,i);
    sh = scatter(ah,GroupVariableNumber(indx),y_var(indx),MarkerSize,ColorMap(i,:),MarkerType{i},'XJitter','density','XJitterWidth',XJitterWidth,'Linewidth',MarkerLineWidth,'MarkerFaceColor',MarkerFaceColor);
    row = dataTipTextRow('',SampleId(indx));
    sh.DataTipTemplate.DataTipRows = row;

end
% Create boxes
bh = boxchart(ah,GroupVariableNumber,y_var,'orientation','vertical','BoxWidth',BoxWidth,'BoxFaceColor','none','BoxEdgeColor',BoxColor,'MarkerStyle','none','LineWidth',BoxLineWidth);

% Add title
if ~isempty(TitleTxt)
    title(ah,TitleTxt,'FontWeight','normal','FontSize',FontSize)
end

% Add Y label
if isempty(YlabelTxt)
    ylabel(sprintf('%s',VariableId{1}),'FontSize',FontSize,'Interpreter','tex')
    %ylabel(sprintf('\\it%s\\rm expression',VariableId{1}),'FontSize',FontSize)
    %ylabel(sprintf('\\it %s\\rm expression',VariableId{1}),'FontSize',FontSize)
    %ylabel(sprintf('\\it%s',VariableId{1}),'FontSize',FontSize)

    %ylabel(sprintf('%s \\beta-value',VariableId{1}),'FontSize',FontSize)
    %ylabel(sprintf('Average \\beta-value',VariableId{1}),'FontSize',FontSize)
else
    ylabel(YlabelTxt,'FontSize',FontSize)
end



% Adjust X-axis
ah.XLim = [0.3 nGroups+0.7];
ah.XTick = 1:nGroups;
ah.XTickLabelRotation = XTickAngle;
ah.XTickLabel = GroupName;

y_nudge=range(y_var)/10;
%ah.YLim = [min(y_var)-y_nudge max(y_var)+y_nudge*2];


% calculate stats and add lines
MAX_Y = max(y_var);
Y_pos = MAX_Y;
if CalcStats
    for i=1:size(CalcGroup,1)
        if 0%CalcGroupAllUnique && diff(CalcGroup(i,:)) == 1
            max_y = max([y_var(SampleIndxMat(:,CalcGroup(i,1))); y_var(SampleIndxMat(:,CalcGroup(i,2)))]);
            Y_pos = max_y+(y_nudge/1.5);
        else
            Y_pos = Y_pos+y_nudge;
        end
        
        %Create line between the two croupd being compared
        line(CalcGroup(i,:),[Y_pos Y_pos],'Color',[0 0 0],'Linewidth',StatLineWidth)
        switch lower(StatType)
            case 't-test'
                [~,p_val] = ttest2(y_var(SampleIndxMat(:,CalcGroup(i,1))),y_var(SampleIndxMat(:,CalcGroup(i,2))),0.05,'both','unequal');
            case 'mw'
                [p_val] = ranksum(y_var(SampleIndxMat(:,CalcGroup(i,1))),y_var(SampleIndxMat(:,CalcGroup(i,2))));
            case 'f-test'
                [~,p_val] = vartest2(y_var(SampleIndxMat(:,CalcGroup(i,1))),y_var(SampleIndxMat(:,CalcGroup(i,2))),0.05);
            case 'bartlett'
                x = y_var(SampleIndxMat(:,CalcGroup(i,1)));
                y = y_var(SampleIndxMat(:,CalcGroup(i,2)));
                p_val = vartestn([x;y],[ones(size(x));ones(size(y))*2],'Display','off','TestType','Bartlett');

        end

        txt_p = pval2stars(p_val,[]);
        if PlotStars
            if p_val > 0.05
                FS = FontSize;
                VerticalAlignment = 'bottom';
            else
                FS = FontSize + 4;
                VerticalAlignment = 'middle';
            end
            if Show_NS || p_val < 0.05
                text(ah,mean(CalcGroup(i,:)),Y_pos+(y_nudge/12),txt_p,'VerticalAlignment',VerticalAlignment,'HorizontalAlignment', 'center','FontSize',FS);
            end
        else
            VerticalAlignment = 'bottom';
            txt_p = num2str(p_val,5);
            text(ah,mean(CalcGroup(i,:)),Y_pos+(y_nudge/12),txt_p,'VerticalAlignment',VerticalAlignment,'HorizontalAlignment', 'center','FontSize',FontSize);
        end


    end
end
try
    ah.YLim = [min(y_var)-y_nudge max([Y_pos,MAX_Y]) + y_nudge*2];
catch

end


