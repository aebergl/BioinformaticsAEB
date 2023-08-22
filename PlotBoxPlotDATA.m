function fh = PlotBoxPlotDATA(DATA,VariableId,GroupVariableName,GroupsToUse,CMap,MarkerTypes,varargin)
MarkerSize = 20;
MarkerLineWidth = 0.5;
BoxLineWidth = 0.5;
FontSize = 10;
FigureSize = [3,2.6];
BoxColor = [0 0 0];
MarkerTypes = {'o'}';
CMap = GetPalette('aeb01');
MarkerFaceColor = 'none';
BoxWidths = 0.8;

VariableIdentifier = false;

CMap = GetPalette('Lancet',[3 4 5]);

SortData=false;
SortData='descend';

i=0;
while i<numel(varargin)
    i = i + 1;
    if strcmpi(varargin{i},'VariableIdentifier')
        i = i + 1;
        VariableIdentifier = varargin{i};
    elseif strcmpi(varargin{i},'Truncate')
        i = i + 1;
        Truncate = varargin{i};
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
indx = ismember(DATA_ID,VariableId);
y_var = DATA.X(:,indx);
y_var(~SampleIndxToUse) = NaN;

% ColorMap
CMap = repmat(CMap,ceil(nGroups/size(CMap,1)),1);
CMap = CMap(1:nGroups,:);
[CMap, descriptorname, description] = colorcet('D8','N',nGroups,'reverse', true);


% Marker Type
if size(MarkerTypes,1) < size(MarkerTypes,2)
    MarkerTypes = MarkerTypes';
end
MarkerTypes = repmat(MarkerTypes,ceil(nGroups/size(MarkerTypes,1)),1);
MarkerTypes = MarkerTypes(1:nGroups,:);


if SortData
    MedianX = zeros(nGroups,1);
    for i=1:nGroups
        indx = SampleIndxMat(:,i);
        MedianX(i) = median(y_var(indx),'omitnan');
    end
    [~,indx] = sort(MedianX,SortData);
    GroupVariableNumberNew = GroupVariableNumber;
    for i=1:nGroups
        GroupVariableNumberNew(GroupVariableNumber == indx(i)) = i;
    end
    GroupVariableNumber= GroupVariableNumberNew;
    %GroupVariableNumber = GroupVariableNumber(indx);
    SampleIndxMat = SampleIndxMat(:,indx);
    GroupName = GroupName(indx);
end

% Create Figure
fh = figure('Name','Box Plot','Color','w','Tag','Box Plot','Units','inches','Colormap',CMap);
fh.Position(3:4) = FigureSize;
ah = axes(fh,'NextPlot','add','tag','Gene Sample Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
    'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');
ah.LineWidth = 0.5;
ah.Colormap=CMap;
for i=1:nGroups
    indx = SampleIndxMat(:,i);
    scatter(ah,GroupVariableNumber(indx),y_var(indx),MarkerSize,CMap(i,:),MarkerTypes{i},'XJitter','density','Linewidth',MarkerLineWidth,'MarkerFaceColor',MarkerFaceColor);
end

bh = boxplot(ah,y_var,GroupVariableNumber,'orientation','vertical','color',BoxColor,'Symbol','','Widths',BoxWidths);
s=findobj( ah.Children, 'Tag', 'boxplot' );
set( findobj( s.Children, 'LineStyle', '--' ),'LineStyle','-')
set( s.Children,'LineWidth',BoxLineWidth)

ylabel(sprintf('\\it %s\\rm expression',VariableId),'FontSize',FontSize)
ah.XLim = [0.3 nGroups+0.7];
ah.XTick = 1:nGroups;
ah.XTickLabelRotation = -45;
ah.XTickLabel = GroupName;

y_nudge=range(y_var)/20;
ah.YLim = [min(y_var)-y_nudge max(y_var)+y_nudge*2];

% line([1 2],[max(y_var)+y_nudge/1.5 max(y_var)+y_nudge/1.5],'Color',[0 0 0],'Linewidth',0.75)
% [~,p_tt] = ttest2(y_var(indx_1),y_var(indx_2),0.05,'both','unequal');
% [p_mw] = ranksum(y_var(indx_1),y_var(indx_2));
% txt_p = pval2stars(p_tt,[]);
% if p_tt > 0.05
%     FontSize = 6;
%     VerticalAlignment = 'baseline';
% else
%     FontSize = 8;
%     VerticalAlignment = 'middle';
% end
%     text(ah,1.5,max(y_var)+y_nudge,txt_p,'VerticalAlignment',VerticalAlignment,'HorizontalAlignment', 'center','FontSize',FontSize);
% 
% 
% UniqueLineObjects=gobjects(nGroups,1);
% for i = 1:nGroups
%     group_indx = find(SampleIndx(:,i));
%     for j = 1:length(group_indx)
%         [N,edges] = histcounts(DATA.X(group_indx(j),:),nBins);
%         x = edges(1:end-1) + diff(edges)/2;
%         UniqueLineObjects(i) = line(x,N,'Parent',ah,'LineWidth',LineWidth, 'LineStyle',MarkerTypes(i),'MarkerEdgeColor',[0 0 0],'Color',CMap(i,:),'MarkerFaceColor',CMap(i,:));
%     end
% end
% if ~isempty(GroupName)
%     lh=legend(UniqueLineObjects,GroupName,'Location','northeastoutside');
%     lh.Interpreter='none';
%     lh.Box='off';
% end
% ah.YTickLabel= [];

