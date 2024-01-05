function fh = PlotScatterDATA(DATA,VariableId_x,VariableId_y,GroupVariableName,GroupsToUse,varargin)
MarkerSize = 100;
MarkerEdgeLineWidth = 0.5;
AlphaValue=0.8;
AlphaValueMarkerLine = 0.8;
MarkerEdgeColor = [0.1 0.1 0.1];
FontSize = 10;
FigureSize = [3,2.6];
MarkerTypes = {'o','v','^','<','>','d'}';
CMap = GetPalette('aeb01');
AxisType = 'normal';
VariableIdentifier = false;

%CMap = GetPalette('Lancet',[1 2 3 4 5]);
MarkerTypes = {'o','^','d'};
i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'VariableIdentifier')
        i = i + 1;
        VariableIdentifier = varargin{i};
    elseif strcmpi(varargin{i},'Truncate')
        i = i + 1;
        Truncate = varargin{i};
    elseif strcmpi(varargin{i},'colormap')
        i = i + 1;
        CMap = varargin{i};
    elseif strcmpi(varargin{i},'MarkerTypes')
        i = i + 1;
        MarkerTypes = varargin{i};

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

% Select Sample Id
SampleId = DATA.RowAnnotation(:,1);
SampleId = DATA.RowId;

% Selection of Y variable
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
CMap = repmat(CMap,ceil(nGroups/size(CMap,1)),1);
CMap = CMap(1:nGroups,:);

% Marker Type
if size(MarkerTypes,1) < size(MarkerTypes,2)
    MarkerTypes = MarkerTypes';
end
MarkerTypes = repmat(MarkerTypes,ceil(nGroups/size(MarkerTypes,1)),1);
MarkerTypes = MarkerTypes(1:nGroups,:);


% Create Figure
fh = figure('Name','Scatter Plot','Color','w','Tag','Scatter Plot','Units','inches','Colormap',CMap);
fh.Position(3:4) = FigureSize;
ah = axes(fh,'NextPlot','add','tag','Gene Sample Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
    'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');
ah.LineWidth = 0.5;
ah.Colormap=CMap;
nudgeX=range(x_var)/20;
ah.XLim = [min(x_var,[],"omitnan")-nudgeX  max(x_var,[],"omitnan") + nudgeX];
nudgeY=range(y_var)/20;
ah.YLim = [min(y_var,[],"omitnan")-nudgeY  max(y_var,[],"omitnan") + nudgeY];

axis(AxisType)

for i=1:nGroups
    indx = SampleIndxMat(:,i);
    sh = scatter(ah,x_var(indx),y_var(indx),MarkerSize,CMap(i,:),MarkerTypes{i},'XJitter','density','Linewidth',MarkerEdgeLineWidth,'MarkerFaceColor','flat','MarkerEdgeColor',MarkerEdgeColor);
    row = dataTipTextRow('',SampleId(indx));
    sh.DataTipTemplate.DataTipRows = row;
    sh.MarkerEdgeAlpha = AlphaValueMarkerLine;
    sh.MarkerFaceAlpha = AlphaValue;

end
if ~isempty(GroupName)
    lh=legend(ah,GroupName,'Location','northeastoutside');
    lh.Interpreter='none';
    lh.Box='off';
    lh.AutoUpdate='off';
end

line(ah,ah.XLim,[0 0],'Linewidth',0.5,'Color','k','LineStyle','-')
line(ah,[0 0],ah.YLim,'Linewidth',0.5,'Color','k','LineStyle','-')
ah.XLabel.String = VariableId_x;
ah.YLabel.String = VariableId_y;




