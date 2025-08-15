function fh = PlotScatterDATA(DATA,VariableId_x,VariableId_y,GroupVariableName,GroupsToUse,varargin)
% USAGE:
%   fh  = PlotScatterDATA(DATA,VariableId,GroupVariableName,GroupsToUse,varargin)
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
%   'VariableIdentifier'   Column Id to be used for selecting variable ['ColId']
%   'FontSize'             FontSize for all text in figure [7]
%   'FigSize'              Vector with figure width and hight in inches [ [2 3.5] ]
%   'AxisType'             Axis Type [ 'square' ]
%   'MarkerSize'           Marker size for scatter points [30]
%   'MarkerLineWidth'      Marker line width for scatter points [1]
%   'MarkerTypes'          Marker types for scatter points,{'o','d','^'} [{'o'}]
%   'AlphaValue'           Marker types for scatter points,{'o','d','^'} [{'o'}]
%   'MarkerLineAlphaValue' Marker types for scatter points,{'o','d','^'} [{'o'}]
%   'MarkerColors'         Marker face colors for scatter points [GetPalette('Science')]
%   'XlabelTxt'            Text to be used for X label, if not defined VariableId will be used
%   'YlabelTxt'            Text to be used for Y label, if not defined VariableId will be used
%   'TitleTxt'             Text for Title
%   'CalcCorr'             Calculate correlation [false]
%   'CorrType'             Correlation Type, Spearman, Pearson, ['Spearman']
%   'CorrLine'             Show correlation line [false]
%   'TargetAxes'           Axes handle for target axist []
%   'DataTipId'            Id to be used for datatip 'RowId' or id from 'RowAnnotationFields' ['RowId']
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% by Anders Berglund, 2025 aebergl at gmail.com                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VariableIdentifier = false; 
FontSize = 12;
FigSize = [5,5];
AxisType = 'square';

MarkerSize = 50;
MarkerLineWidth = 0.2;
MarkerTypes = {'o','v','^','<','>','d'}';
AlphaValue = 0.5;
MarkerLineAlphaValue = 0.95;
MarkerColors = GetPalette('aeb01');
MarkerEdgeColor = [0 0 0];
XlabelTxt = [];
YlabelTxt = [];
TitleTxt = [];

CalcCorr = true;
CorrType = 'Spearman';
CorrLine = false;
CorrLineWidth = 1;
CorrLineColor = [0 0 0];
CorrLineStyle = '-';

TargetAxes = false;
DataTipId = 'RowId';


% Check input
if nargin > 5
    ArgsList = {'VariableIdentifier','FontSize','FigSize','AxisType',...
        'MarkerSize','MarkerLineWidth','MarkerTypes','AlphaValue','MarkerLineAlphaValue',...
        'MarkerColors','MarkerEdgeColor','YlabelTxt','YlabelTxt','TitleText',...
        'CalcCorr','CorrType','CorrLine','CorrLineWidth','CorrLineColor','CorrLineType',...
        'TargetAxes','DataTipId'};
    for j=1:2:numel(varargin)
        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)
                case 'variableidentifier'
                    VariableIdentifier = ArgVal;
                 case 'fontsize'
                    FontSize = ArgVal;
                case 'figsize'
                    FigSize = ArgVal;
                case 'axistype'
                    AxisType = ArgVal;
                case 'markersize'
                    MarkerSize = ArgVal;
                case 'markerlinewidth'
                    MarkerLineWidth = ArgVal;
                case 'markertypes'
                    MarkerTypes = ArgVal;
                case 'alphavalue'
                    AlphaValue = ArgVal;
                case 'markerlinealphavalue'
                    MarkerLineAlphaValue = ArgVal;
                case 'markeredgecolor'
                    MarkerEdgeColor = ArgVal;
                case 'markercolors'
                    MarkerColors = ArgVal;
                case 'xlabel'
                    XlabelTxt =ArgVal;
                case 'ylabel'
                    YlabelTxt =ArgVal;
                case 'titletext'
                    TitleTxt = ArgVal;
                case 'calccorr'
                    CalcCorr = ArgVal;
                case 'corrtype'
                    CorrType = ArgVal;
                case 'corrline'
                    CorrLine = true;
                case 'corrlinecolor'
                    CorrLineColor = ArgVal;
                case 'corrlinestyle'
                    CorrLineStyle = ArgVal;
                case 'targetaxes'
                    TargetAxes = ArgVal;
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
SampleIndxToUse = any(SampleIndxMat,2);

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

% Selection of X Y variable id
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

indx = ismember(DATA_ID,VariableId_x);
x_var = DATA.X(:,indx);
x_var(~SampleIndxToUse) = NaN;

indx = ismember(DATA_ID,VariableId_y);
y_var = DATA.X(:,indx);
y_var(~SampleIndxToUse) = NaN;

% ColorMap
MarkerColors = repmat(MarkerColors,ceil(nGroups/size(MarkerColors,1)),1);
MarkerColors = MarkerColors(1:nGroups,:);

% Marker Type
if size(MarkerTypes,1) < size(MarkerTypes,2)
    MarkerTypes = MarkerTypes';
end
MarkerTypes = repmat(MarkerTypes,ceil(nGroups/size(MarkerTypes,1)),1);
MarkerTypes = MarkerTypes(1:nGroups,:);


% Create Figure

% Create Figure
if isgraphics(TargetAxes,'axes')
    ah = TargetAxes;
    fh = TargetAxes.Parent;
else
    fh = figure('Name','Box Plot','Color','w','Tag','Box Plot','Units','inches','Colormap',MarkerColors);
    fh.Position(3:4) = FigSize;
    fh.Renderer='painters';
    ah = axes(fh,'NextPlot','add','tag','Gene Sample Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
        'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');
    ah.LineWidth = 0.5;
end
ah.Colormap=MarkerColors;

nudgeX=range(x_var)/20;
ah.XLim = [min(x_var,[],"omitnan")-nudgeX  max(x_var,[],"omitnan") + nudgeX];
nudgeY=range(y_var)/20;
ah.YLim = [min(y_var,[],"omitnan")-nudgeY  max(y_var,[],"omitnan") + nudgeY];

axis(AxisType)

for i=1:nGroups
    indx = SampleIndxMat(:,i);
    sh(i) = scatter(ah,x_var(indx),y_var(indx),MarkerSize,MarkerColors(i,:),MarkerTypes{i},'MarkerFaceColor','flat','MarkerEdgeColor',MarkerEdgeColor);
    row = dataTipTextRow('',SampleId(indx));
    sh(i).DataTipTemplate.DataTipRows = row;
    sh(i).DataTipTemplate.Interpreter='none';
    sh(i).MarkerEdgeAlpha = MarkerLineAlphaValue;
    sh(i).MarkerFaceAlpha = AlphaValue;
        sh(i).LineWidth = MarkerLineWidth;


end
if ~isempty(GroupName)
    lh=legend(ah,GroupName,'Location','northeastoutside');
    lh.Interpreter='none';
    lh.Box='off';
    lh.AutoUpdate='off';
end

if CalcCorr
    [r, p_val] = corr(x_var, y_var,'Type',CorrType,'Rows','pairwise');
    txt_str{1} = sprintf('r=%.3f',r);
    txt_str{2} = sprintf('p=%.3g',p_val);
    text(ah,(ah.XLim(1))+nudgeX,ah.YLim(2)-nudgeY/2,txt_str,'VerticalAlignment','Top','HorizontalAlignment', 'Left','FontSize',FontSize);
end
if CorrLine
    p = polyfit(x_var,y_var,1);
    f = polyval(p,x_var);
%    hold on
    line(ah,x_var,f,'Color',CorrLineColor,'LineStyle',CorrLineStyle,'LineWidth',CorrLineWidth);
end


% Add title
if ~isempty(TitleTxt)
    title(ah,TitleTxt,'FontWeight','normal','FontSize',FontSize)
end

% Add Y label
if isempty(XlabelTxt)
    %xlabel(sprintf('%s',VariableId_x),'FontSize',FontSize,'Interpreter','tex')
    xlabel(sprintf('\\it%s\\rm expression',VariableId_x),'FontSize',FontSize)
else
    xlabel(YlabelTxt,'FontSize',FontSize)
end

if isempty(YlabelTxt)
    ylabel(sprintf('%s',VariableId_y),'FontSize',FontSize,'Interpreter','tex')
    %ylabel(sprintf('\\it%s\\rm expression',VariableId_y),'FontSize',FontSize)

else
    ylabel(YlabelTxt,'FontSize',FontSize)
end




