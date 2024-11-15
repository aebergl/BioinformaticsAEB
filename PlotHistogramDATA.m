function fh = PlotHistogramDATA(DATA,GroupVariable,GroupsToUse,CMap,LineTypes)
LineWidth = 0.5;
FontSize = 12;
nBins = 100;
BandwidthValue = 0.05;
nPoints=1000;
KernalDensity = false;
AlphaValue  = 0.8;

if isempty(GroupVariable)
    nGroups = 1;
    SampleIndx = true(DATA.nRow,nGroups);
    GroupName = [];
else
    indx = strcmpi(GroupVariable,DATA.RowAnnotationFields);
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
    SampleIndx = true(DATA.nRow,nGroups);
    GroupVariable = cell(DATA.nRow,1);
    GroupName = cell(nGroups,1);
    for i = 1:nGroups
        tmp_name = GroupsToUse{i};
        if ~iscell(tmp_name)
            tmp_name = {tmp_name};
        end
        SampleIndxMatrix = false(DATA.nRow,numel(tmp_name));
        for j=1:length(tmp_name)
            SampleIndxMatrix(:,j) = strcmp(tmp_name{j},GroupData);
        end
        SampleIndx(:,i) = any(SampleIndxMatrix,2);
        GroupVariable(SampleIndx(:,i)) = tmp_name(1);
        numValuesPerGroup(i) = sum(SampleIndx(:,i));
        GroupName(i) = tmp_name(1);
    end
end

if isempty(CMap)
    CMap = GetPalette('Tab10');
end
CMap = repmat(CMap,ceil(nGroups/size(CMap,1)),1);
CMap = CMap(1:nGroups,:);

if size(LineTypes,1) < size(LineTypes,2)
    LineTypes = LineTypes';
end
if isempty(LineTypes)
    LineTypes = {'-'}';
end
LineTypes = repmat(LineTypes,ceil(nGroups/size(LineTypes,1)),1);
LineTypes = LineTypes(1:nGroups,:);

fh = figure('Name','Histogram Plot','Color','w','Tag','Histogram Plot','Units','inches','Colormap',CMap);

ah = axes(fh,'NextPlot','add','tag','Histogram Plot','Box','on','FontSize',FontSize,'Linewidth',0.5,...
    'ActivePositionProperty','outerposition','XGrid','on','YGrid','on');

UniqueLineObjects=gobjects(nGroups,1);

for i = 1:nGroups
    group_indx = find(SampleIndx(:,i));
    for j = 1:length(group_indx)
        if KernalDensity 
            drawnow
            [N,x] = ksdensity(DATA.X(group_indx(j),:),'npoints',nPoints,'Bandwidth',BandwidthValue);

        else
            [N,edges] = histcounts(DATA.X(group_indx(j),:),nBins);
            x = edges(1:end-1) + diff(edges)/2;


        end

        UniqueLineObjects(i) = line(x,N,'Parent',ah,'LineWidth',LineWidth,'Color',[CMap(i,:) AlphaValue],'LineStyle',LineTypes{i});
        
        % row = dataTipTextRow('',DATA.RowAnnotation(group_indx(j),1));
        % UniqueLineObjects(i).DataTipTemplate.DataTipRows = row;
        % UniqueLineObjects(i).DataTipTemplate.Interpreter='none';
    end
end
if ~isempty(GroupName)
    lh=legend(UniqueLineObjects,GroupName,'Location','northeastoutside');
    lh.Interpreter='none';
    lh.Box='off';
end
ah.YTickLabel= [];
%ah.XLim=[-0.05 1.05];
% ah.XTick = 1:nProbes;
% ah.XTickLabelRotation=-45;
% text_str = DATA.ColAnnotation(ProbeIndx,2);
% text_str = strcat('\it',text_str);
% ah.XTickLabel = text_str;
% 
% ylabel('Log_2 Expression Value')
