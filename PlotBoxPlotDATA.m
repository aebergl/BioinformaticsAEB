function fh = PlotBoxPlotDATA(DATA,VariableId,GroupVariableName,GroupsToUse,varargin)
MarkerSize = 10;
MarkerLineWidth = 0.5;
BoxLineWidth = 0.5;
FontSize = 6;
FigureSize = [3,2.6];
BoxColor = [0 0 0];
MarkerTypes = {'o'}';
CMap = GetPalette('aeb01');
MarkerFaceColor = 'none';
BoxWidths = 0.8;
XJitterWidth = 0.6;
YlabelTxt = [];
TitleTxt = [];
VariableIdentifier = false;
CalcStats = false;
CalcGroup = [];
CalcGroupAllUnique = true;
StatType = 't-test';
StatType = 'MW';
PlotStars = true;
TargetAxes = false;
XTickAngle = -45;
% CMap = GetPalette('Lancet',[3 4 5]);

CMap = GetPalette('Science');

SortData=false;
%SortData='descend';

% Check input
if nargin > 4
    ArgsList = {'VariableIdentifier','FigureSize','Ylabel','TitleText','ColorMap','MarkerTypes','CalcStats','TargetAxes'};
    for j=1:2:numel(varargin)
        ArgType = varargin{j};
        ArgVal = varargin{j+1};
        if ~strncmpi(ArgType,ArgsList,numel(ArgType))
            error('Invalid input option: %s', ArgType);
        else
            switch lower(ArgType)

                case 'variableidentifier'
                    VariableIdentifier =ArgVal;
                case 'figuresize'
                    FigureSize = ArgVal;
                case 'ylabel'
                    YlabelTxt =ArgVal;
                case 'titletext'
                    TitleTxt = ArgVal;
                case 'colormap'
                    CMap = ArgVal;
                case 'markertypes'
                    MarkerTypes = ArgVal;
                case 'calcstats'
                    CalcGroup = ArgVal;
                    CalcStats = true;
               case 'targetaxes'
                    TargetAxes = ArgVal;

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

if size(y_var,2) > 1
    %y_var = B2M(y_var);
    %y_var=mean(y_var,2,"omitnan");
    [PCA_Model] = NIPALS_PCA(y_var,'NumComp',1,'ScaleX',false);
    y_var = PCA_Model.T(:,1);
    VariableId = "Average";
end

% ColorMap
CMap = repmat(CMap,ceil(nGroups/size(CMap,1)),1);
CMap = CMap(1:nGroups,:);
%[CMap, descriptorname, description] = colorcet('D8','N',nGroups,'reverse', true);


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
if isgraphics(TargetAxes,'axes')
    ah = TargetAxes;
    fh = TargetAxes.Parent;
else
    fh = figure('Name','Box Plot','Color','w','Tag','Box Plot','Units','inches','Colormap',CMap);
    fh.Position(3:4) = FigureSize;
    ah = axes(fh,'NextPlot','add','tag','Gene Sample Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
        'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');
    ah.LineWidth = 0.5;
end
ah.Colormap=CMap;

% Select Sample Id
%SampleId = DATA.RowAnnotation(:,18);
SampleId = DATA.RowId;
for i=1:nGroups
    indx = SampleIndxMat(:,i);
    sh = scatter(ah,GroupVariableNumber(indx),y_var(indx),MarkerSize,CMap(i,:),MarkerTypes{i},'XJitter','density','XJitterWidth',XJitterWidth,'Linewidth',MarkerLineWidth,'MarkerFaceColor',MarkerFaceColor);
    row = dataTipTextRow('',SampleId(indx));
    sh.DataTipTemplate.DataTipRows = row;

end

%bh = boxplot(ah,y_var,GroupVariableNumber,'orientation','vertical','color',BoxColor,'Symbol','','Widths',BoxWidths);
bh = boxchart(ah,GroupVariableNumber,y_var,'orientation','vertical','BoxWidth',BoxWidths,'BoxFaceColor','none','BoxEdgeColor',BoxColor,'MarkerStyle','none','LineWidth',BoxLineWidth);

%s=findobj( ah.Children, 'Tag', 'boxchart' )
% set( findobj( s.Children, 'LineStyle', '--' ),'LineStyle','-')
% set( s.Children,'LineWidth',BoxLineWidth)

if isempty(YlabelTxt)
    %ylabel(sprintf('\\it %s\\rm value',VariableId{1}),'FontSize',FontSize)
    ylabel(sprintf('\\it%s',VariableId{1}),'FontSize',FontSize)
    %ylabel(sprintf('%s \\beta-value',VariableId{1}),'FontSize',FontSize)
    %ylabel(sprintf('Average \\beta-value',VariableId{1}),'FontSize',FontSize)
else
    ylabel(YlabelTxt,'FontSize',FontSize)
end

if ~isempty(TitleTxt)
    title(ah,TitleTxt,'FontWeight','normal','FontSize',FontSize)
end

ah.XLim = [0.3 nGroups+0.7];
ah.XTick = 1:nGroups;
ah.XTickLabelRotation = XTickAngle;
ah.XTickLabel = GroupName;

y_nudge=range(y_var)/10;
%ah.YLim = [min(y_var)-y_nudge max(y_var)+y_nudge*2];


MAX_Y = max(y_var);
Y_pos = MAX_Y;
if CalcStats
    for i=1:size(CalcGroup,1)
        if CalcGroupAllUnique && diff(CalcGroup(i,:)) == 1
            max_y = max([y_var(SampleIndxMat(:,CalcGroup(i,1))); y_var(SampleIndxMat(:,CalcGroup(i,2)))]);
            Y_pos = max_y+y_nudge/1.5;
        else
            Y_pos = Y_pos+y_nudge;
        end
        line(CalcGroup(i,:),[Y_pos Y_pos],'Color',[0 0 0],'Linewidth',0.75)
        switch lower(StatType)
            case 't-test'
                [~,p_val] = ttest2(y_var(SampleIndxMat(:,CalcGroup(i,1))),y_var(SampleIndxMat(:,CalcGroup(i,2))),0.05,'both','unequal');
            case 'mw'
                [p_val] = ranksum(y_var(SampleIndxMat(:,CalcGroup(i,1))),y_var(SampleIndxMat(:,CalcGroup(i,2))));
        end
        
        txt_p = pval2stars(p_val,[]);
        if PlotStars
            if p_val > 0.05
                FontSize = FontSize;
                VerticalAlignment = 'bottom';
            else
                FontSize = 8;
                VerticalAlignment = 'middle';
            end
            text(ah,mean(CalcGroup(i,:)),Y_pos+(y_nudge/12),txt_p,'VerticalAlignment',VerticalAlignment,'HorizontalAlignment', 'center','FontSize',FontSize);
        else

            FontSize = 6;
            VerticalAlignment = 'bottom';
            txt_p = num2str(p_val,2);
            text(ah,mean(CalcGroup(i,:)),Y_pos+(y_nudge/12),txt_p,'VerticalAlignment',VerticalAlignment,'HorizontalAlignment', 'center','FontSize',FontSize);
        end


    end
end
ah.YLim = [min(y_var)-y_nudge max([Y_pos,MAX_Y]) + y_nudge*2];
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

